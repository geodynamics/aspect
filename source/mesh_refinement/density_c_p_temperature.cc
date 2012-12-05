/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
/*  $Id$  */


#include <aspect/mesh_refinement/density_c_p_temperature.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/derivative_approximation.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    DensityCpTemperature<dim>::execute(Vector<float> &indicators) const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("This refinement criterion can not be used when no "
                               "compositional fields are active!"));
      indicators = 0;

      // create a vector in which we set the temperature block to
      // be a finite element interpolation of the density or rho*c_p*T.
      // we do so by setting up a quadrature formula with the
      // temperature unit support points, then looping over these
      // points, compute the output quantity at them, and writing
      // the result into the output vector in the same order
      // (because quadrature points and temperature dofs are,
      // by design of the quadrature formula, numbered in the
      // same way)
//TODO: There should be a way to get at this kind of information via SimulatorAccess
      std::vector<unsigned int> system_sub_blocks (dim+2+this->n_compositional_fields(),0);
      system_sub_blocks[dim] = 1;
      system_sub_blocks[dim+1] = 2;
      for (unsigned int i=dim+2; i<dim+2+this->n_compositional_fields(); ++i)
        system_sub_blocks[i] = i-dim+1;
      std::vector<unsigned int> system_dofs_per_block (3+this->n_compositional_fields());
      DoFTools::count_dofs_per_block (this->get_dof_handler(), system_dofs_per_block,
                                      system_sub_blocks);

      const unsigned int n_u = system_dofs_per_block[0],
                         n_p = system_dofs_per_block[1],
                         n_T = system_dofs_per_block[2];
      unsigned int       n_C_sum = 0;
      std::vector<unsigned int> n_C (this->n_compositional_fields()+1);
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          n_C[c] = system_dofs_per_block[c+3];
          n_C_sum += n_C[c];
        }

      std::vector<IndexSet> system_partitioning;
      {
        IndexSet system_index_set = this->get_dof_handler().locally_owned_dofs();
        system_partitioning.push_back(system_index_set.get_view(0,n_u));
        system_partitioning.push_back(system_index_set.get_view(n_u,n_u+n_p));
        system_partitioning.push_back(system_index_set.get_view(n_u+n_p,n_u+n_p+n_T));

        {
          unsigned int n_C_so_far = 0;

          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              system_partitioning.push_back(system_index_set.get_view(n_u+n_p+n_T+n_C_so_far,
                                                                      n_u+n_p+n_T+n_C_so_far+n_C[c]));
              n_C_so_far += n_C[c];
            }
        }
      }
      LinearAlgebra::BlockVector vec_distributed (system_partitioning, this->get_mpi_communicator());

      const Quadrature<dim> quadrature(this->get_fe().base_element(2).get_unit_support_points());
      std::vector<unsigned int> local_dof_indices (this->get_fe().dofs_per_cell);
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values);
      std::vector<double> pressure_values(quadrature.size());
      std::vector<double> temperature_values(quadrature.size());

      // the values of the compositional fields are stored as blockvectors for each field
      // we have to extract them in this structure
      std::vector<std::vector<double> > prelim_composition_values (this->n_compositional_fields(),
                                                                   std::vector<double> (quadrature.size()));
      std::vector<std::vector<double> > composition_values (quadrature.size(),
                                                            std::vector<double> (this->n_compositional_fields()));

      const FEValuesExtractors::Scalar pressure (dim);
      const FEValuesExtractors::Scalar temperature (dim+1);
      std::vector<FEValuesExtractors::Scalar> composition;

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          const FEValuesExtractors::Scalar temp(dim+2+c);
          composition.push_back(temp);
        }

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values[pressure].get_function_values (this->get_solution(),
                                                     pressure_values);
            fe_values[temperature].get_function_values (this->get_solution(),
                                                        temperature_values);
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[composition[c]].get_function_values (this->get_solution(),
                                                             prelim_composition_values[c]);
            // then we copy these values to exchange the inner and outer vector, because for the material
            // model we need a vector with values of all the compositional fields for every quadrature point
            for (unsigned int q=0; q<quadrature.size(); ++q)
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                composition_values[q][c] = prelim_composition_values[c][q];

            cell->get_dof_indices (local_dof_indices);

            // for each temperature dof, write into the output
            // vector the density. note that quadrature points and
            // dofs are enumerated in the same order
            for (unsigned int i=0; i<this->get_fe().base_element(2).dofs_per_cell; ++i)
              {
                const unsigned int system_local_dof
                  = this->get_fe().component_to_system_index(/*temperature component=*/dim+1,
                                                                                       /*dof index within component=*/i);

                vec_distributed(local_dof_indices[system_local_dof])
                  = this->get_material_model().density( temperature_values[i],
                                                        pressure_values[i],
                                                        composition_values[i],
                                                        fe_values.quadrature_point(i))
                    * temperature_values[i]
                    * this->get_material_model().specific_heat(temperature_values[i],
                                                               pressure_values[i],
                                                               composition_values[i],
                                                               fe_values.quadrature_point(i));
              }
          }

      // now create a vector with the requisite ghost elements
      // and use it for estimating the gradients
      LinearAlgebra::BlockVector vec (this->get_solution());
      vec = vec_distributed;

      DerivativeApproximation::approximate_gradient  (this->get_mapping(),
                                                      this->get_dof_handler(),
                                                      vec,
                                                      indicators,
                                                      dim+1);

      // Scale gradient in each cell with the correct power of h. Otherwise,
      // error indicators do not reduce when refined if there is a density
      // jump. We need at least order 1 for the error not to grow when
      // refining, so anything >1 should work.
      const double power = 1.5;
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();
        unsigned int i=0;
        for (; cell!=endc; ++cell, ++i)
          if (cell->is_locally_owned())
            indicators(i) *= std::pow(cell->diameter(), power);
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(DensityCpTemperature,
                                              "density cp temperature",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators from a field that describes "
                                              "the spatial variability of the thermal energy density, $\\rho C_p T$. "
                                              "Unlike most other criteria, this one uses the gradient, rather than "
                                              "the second derivative of the energy density field to compute error "
                                              "indicators.")
  }
}
