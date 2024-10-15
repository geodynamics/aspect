/*
  Copyright (C) 2020 - 2024 by the authors of the ASPECT code.

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

#include <aspect/postprocess/visualization/boundary_strain_rate_residual.h>
#include <aspect/postprocess/boundary_strain_rate_residual_statistics.h>
#include <aspect/simulator.h>

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
        DataPostprocessorScalar<dim> ("boundary_strain_rate_residual",
                                      update_quadrature_points | update_gradients),
        Interface<dim>()    // unknown units at construction time, will be filled by a separate function
      {}



      template <int dim>
      std::string
      BoundaryStrainRateResidual<dim>::
      get_physical_units () const
      {
        if (this->convert_output_to_years())
          return "1/year";
        else
          return "1/s";
      }



      template <int dim>
      void
      BoundaryStrainRateResidual<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        (void) n_quadrature_points;
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError());
        Assert (computed_quantities[0].size() == 1,               ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components, ExcInternalError());

        auto cell = input_data.template get_cell<dim>();

        for (auto &quantity : computed_quantities)
          quantity(0)= std::numeric_limits<double>::quiet_NaN();

        const double unit_scaling_factor = this->convert_output_to_years() ? year_in_seconds : 1.0;

        const Postprocess::BoundaryStrainRateResidualStatistics<dim> &boundary_strain_rate_residual_statistics =
          this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::BoundaryStrainRateResidualStatistics<dim>>();

        // We only want the output at the top boundary, so only compute it if the current cell
        // has a face at the top boundary.
        bool cell_at_top_boundary = false;
        for (const unsigned int f : cell->face_indices())
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            {
              cell_at_top_boundary = true;
              break;
            }

        if (cell_at_top_boundary)
          for (unsigned int q=0; q<computed_quantities.size(); ++q)
            {
              Tensor<2,dim> grad_u;
              for (unsigned int d=0; d<dim; ++d)
                grad_u[d] = input_data.solution_gradients[q][d];
              const SymmetricTensor<2,dim> strain_rate = symmetrize(grad_u);

              const double data_surface_strain_rate =
                boundary_strain_rate_residual_statistics.get_data_surface_strain_rate(input_data.evaluation_points[q]);

              // Only compute residual for doubles. This condition checks for nan values
              if (data_surface_strain_rate < 1e300 || !std::isnan(data_surface_strain_rate))
                computed_quantities[q](0) = data_surface_strain_rate -
                                            std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))) * unit_scaling_factor;
              else
                continue;

            }

        const auto &viz = this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Visualization<dim>>();
        if (!viz.output_pointwise_stress_and_strain())
          average_quantities(computed_quantities);

      }



      template <int dim>
      std::list<std::string>
      BoundaryStrainRateResidual<dim>::required_other_postprocessors() const
      {
        return {"boundary strain rate residual statistics"};
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
                                                  "A visualization output object that generates output for the strain rate "
                                                  "residual at the top surface. The residual is computed at each point at the "
                                                  "surface as the difference between the strain rate invariant in the model and the input data, "
                                                  "where the invariant is computed like in the 'strain rate' postprocessor. The user chooses "
                                                  "the input data as ascii data files with coordinate columns and column corresponding "
                                                  "to the surface strain rate norm."
                                                  "\n\n"
                                                  "Physical units: $\\frac{1}{\\text{s}}$ or "
                                                  "$\\frac{1}{\\text{year}}$, depending on settings in the input file.")
    }
  }
}
