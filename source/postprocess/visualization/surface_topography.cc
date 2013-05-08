/*
  Copyright (C) 2012 by the authors of the ASPECT code.

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
/*  $Id: surface_topography.cc 1298 2012-10-19 15:41:38Z dannberg $  */


#include <aspect/postprocess/visualization/surface_topography.h>
#include <aspect/simulator_access.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
    namespace
    {
      /**
       * Copy the values of the compositional fields at the quadrature point
       * q given as input parameter to the output vector
       * composition_values_at_q_point.
       */
      void
      extract_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                             const unsigned int                      q,
                                             std::vector<double>                    &composition_values_at_q_point)
      {
        for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
          composition_values_at_q_point[k] = composition_values[k][q];
      }
    }

      template <int dim>
      std::pair<std::string, Vector<float> *>
      SurfaceTopography<dim>::execute() const
      {
    	// check that the topography refinement criterion is switched on
    	// otherwise the result will be inaccurate, because different cell
    	// sizes lead to the evaluation of the pressure in different depths

        std::pair<std::string, Vector<float> *>
        return_value ("topography",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        // evaluate a single point per cell
        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_gradients   |
                                 update_quadrature_points );

        std::vector<double> pressure_values(n_q_points);
        std::vector<double> temperature_values(n_q_points);
        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),
                                                              std::vector<double> (n_q_points));
        std::vector<double> composition_values_at_q_point (this->n_compositional_fields());
        std::vector<SymmetricTensor<2,dim> > strain_rate(n_q_points);

        double local_topography_integral = 0;
        double local_length_integral = 0;

        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              // get all the variables we need
              fe_values.reinit (cell);
              fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                        pressure_values);
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                           temperature_values);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    composition_values[c]);
              fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients(this->get_solution(),
                                                                                                      strain_rate);

              extract_composition_values_at_q_point (composition_values,
                                                     0,
                                                     composition_values_at_q_point);

              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));

              // construct normal vector
              Tensor<1,dim> normal(0.0);
              if (dynamic_cast <const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0)
                for (unsigned int i=0;i<dim;++i)
            	  normal[i] = fe_values.quadrature_point(0)[i]/fe_values.quadrature_point(0).norm();
              else if (dynamic_cast <const GeometryModel::Box<dim>*> (&this->get_geometry_model()) != 0)
                normal[dim-1] = 1.0;
              else ExcMessage ("Not a valid geometry model for the postprocessor"
            		             "surface topography.");

              double cell_depth = 0;
              for (unsigned int i=0;i<dim;++i)
            	cell_depth += std::abs(normal[i])*cell->extent_in_direction(i);

              // only calculate topography for upper surface
              if (cell->at_boundary() && depth < cell_depth)
              {
                // evaluate material properties
                typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, this->n_compositional_fields());
                typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1);
                in.position[0] = fe_values.quadrature_point(0);
                in.temperature[0] = temperature_values[0];
                in.pressure[0] = pressure_values[0];
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[0][c] = composition_values_at_q_point[c];
                in.strain_rate[0] = strain_rate[0];

                this->get_material_model().evaluate(in, out);

                const SymmetricTensor<2,dim> compressible_strain_rate
                  = (this->get_material_model().is_compressible()
                     ?
                     strain_rate[0] - 1./3 * trace(strain_rate[0]) * unit_symmetric_tensor<dim>()
                     :
                     strain_rate[0]);

                double shear_stress = 0.0;
                const double viscosity = out.viscosities[0];

                Tensor<2,dim> strain_rate_tensor;
                for (unsigned int i=0;i<dim;++i)
                  for (unsigned int j=0;j<dim;++j)
                    strain_rate_tensor[i][j] = compressible_strain_rate[i][j];

                // get the normal stress
                Tensor<1,dim> normal_strain_rate;
                contract(normal_strain_rate,strain_rate_tensor,normal);

                shear_stress = 2.0 * viscosity * normal_strain_rate * normal;

                const double dynamic_pressure = pressure_values[0] - this->get_adiabatic_conditions().pressure(fe_values.quadrature_point(0));
                const double normal_stress = - dynamic_pressure + shear_stress;

                const double density = out.densities[0];
                const double surface_density = 0.0;
                const double gravity = this->get_gravity_model().gravity_vector(fe_values.quadrature_point(0)).norm();

                // also calculate average topography for normalization
                (*return_value.second)(cell_index) = - normal_stress / ((density - surface_density) * gravity);
                local_topography_integral += (*return_value.second)(cell_index)*cell->diameter();
                local_length_integral += cell->diameter();
              }
              else
                (*return_value.second)(cell_index) = 0;
            }
        // calculate global average over all processors
        const double global_length_integral
          = Utilities::MPI::sum (local_length_integral, this->get_mpi_communicator());
        const double global_topography_average
          = Utilities::MPI::sum (local_topography_integral, this->get_mpi_communicator())/global_length_integral;

        // subtract global average
        cell_index = 0;
        cell = this->get_dof_handler().begin_active();
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned() && (*return_value.second)(cell_index) != 0)
            if (cell->at_boundary())
              (*return_value.second)(cell_index) -= global_topography_average;

        return return_value;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceTopography,
                                                  "surface topography",
												  "A visualization output object that generates output "
												  "for the surface topography using the normal stresses"
												  "at the surface for the calculation. ")
    }
  }
}
