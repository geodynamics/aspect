/*
  Copyright (C) 2020 - 2021 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/simulator.h>
#include <aspect/postprocess/boundary_velocity_residual_statistics.h>
#include <aspect/postprocess/visualization/boundary_strain_rate_residual.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      BoundaryStrainRateResidual<dim>::
      BoundaryStrainRateResidual ()
        :
        DataPostprocessorTensor<dim> ("boundary_strain_rate_residual",
                                      update_quadrature_points | update_gradients)
      {}



      template <int dim>
      void
      BoundaryStrainRateResidual<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert ((computed_quantities[0].size() == Tensor<2,dim>::n_independent_components),
                ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

#if DEAL_II_VERSION_GTE(9,3,0)
        auto cell = input_data.template get_cell<dim>();
#else
        auto cell = input_data.template get_cell<DoFHandler<dim> >();
#endif

        const double scaling_factor =
          this->convert_output_to_years() ? year_in_seconds : 1.0;


        const Postprocess::BoundaryVelocityResidualStatistics<dim> &boundary_velocity_residual_statistics =
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::BoundaryVelocityResidualStatistics<dim> >();

        // We only want the output at the top boundary, so only compute it if the current cell
        // has a face at the top boundary.
        bool cell_at_top_boundary = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            cell_at_top_boundary = true;

        if (cell_at_top_boundary)
          for (unsigned int q=0; q<computed_quantities.size(); ++q)
            {
              Tensor<2,dim> grad_u;
              const std::vector<Vector<double> > data_surface_strain_rate =
                boundary_velocity_residual_statistics.get_data_surface_strain_rate(input_data.evaluation_points[q]);

              // Get the surface coordinates only from the computed values
              if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::cartesian)
                {
                  grad_u[0] = input_data.solution_gradients[q][0];
                  grad_u[1] = input_data.solution_gradients[q][1];
                  grad_u[2] = input_data.solution_gradients[q][2];
                }
              else if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
                {
                  grad_u[0] = input_data.solution_gradients[q][1];
                  grad_u[1] = input_data.solution_gradients[q][2];
                  grad_u[2] = input_data.solution_gradients[q][0];
                }

              const SymmetricTensor<2,dim> strain_rate = symmetrize(grad_u);

              for (unsigned int d=0; d<dim-1; ++d)
                for (unsigned int e=0; e<dim-1; ++e)
                  computed_quantities[q][Tensor<2,dim-1>::component_to_unrolled_index(TableIndices<2>(d,e))] =
                    data_surface_strain_rate[d][e] - strain_rate[d][e] * scaling_factor;

            }

      }



      template <int dim>
      std::list<std::string>
      BoundaryStrainRateResidual<dim>::required_other_postprocessors() const
      {
        return std::list<std::string> (1, "boundary velocity residual statistics");
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(BoundaryStrainRateResidual,
                                                  "boundary strain rate residual",
                                                  "A visualization output object that generates output for the velocity "
                                                  "residual at the top surface. The residual is computed at each point at the "
                                                  "surface as the difference between the modeled velocities and the input "
                                                  "data velocities for each vector component. The user has an option to choose "
                                                  "the input data as ascii data files (e.g. GPS velocities) with columns "
                                                  "in the same format as described for the 'ascii data' initial temperature plugin "
                                                  "or a velocity field computed from the GPlates program as described in the gplates "
                                                  "boundary velocity plugin. ")
    }
  }
}
