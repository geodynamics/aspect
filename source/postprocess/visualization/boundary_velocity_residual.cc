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

#include <aspect/postprocess/visualization/boundary_velocity_residual.h>
#include <aspect/postprocess/boundary_velocity_residual_statistics.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      BoundaryVelocityResidual<dim>::
      BoundaryVelocityResidual ()
        :
        DataPostprocessorVector<dim> ("boundary_velocity_residual",
                                      update_values | update_quadrature_points | update_gradients),
        Interface<dim>()  // unknown units at construction time, will be filled by a separate function
      {}



      template <int dim>
      std::string
      BoundaryVelocityResidual<dim>::
      get_physical_units () const
      {
        if (this->convert_output_to_years())
          return "m/year";
        else
          return "m/s";
      }



      template <int dim>
      void
      BoundaryVelocityResidual<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        Assert ((computed_quantities[0].size() == dim), ExcInternalError());
        for (auto &quantity : computed_quantities)
          for (unsigned int d = 0; d < dim; ++d)
            quantity(d)= 0.;

        auto cell = input_data.template get_cell<dim>();
        const double velocity_scaling_factor =
          this->convert_output_to_years() ? year_in_seconds : 1.0;

        const Postprocess::BoundaryVelocityResidualStatistics<dim> &boundary_velocity_residual_statistics =
          this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::BoundaryVelocityResidualStatistics<dim>>();

        // We only want the output at the top boundary, so only compute it if the current cell
        // has a face at the top boundary.
        bool cell_at_top_boundary = false;
        for (const unsigned int f : cell->face_indices())
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            cell_at_top_boundary = true;

        if (cell_at_top_boundary)
          for (unsigned int q=0; q<computed_quantities.size(); ++q)
            {
              const Tensor<1,dim> data_velocity = boundary_velocity_residual_statistics.get_data_velocity(input_data.evaluation_points[q]);
              for (unsigned int d = 0; d < dim; ++d)
                computed_quantities[q](d) = data_velocity[d] - input_data.solution_values[q][d] * velocity_scaling_factor;
            }

      }



      template <int dim>
      std::list<std::string>
      BoundaryVelocityResidual<dim>::required_other_postprocessors() const
      {
        return {"boundary velocity residual statistics"};
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(BoundaryVelocityResidual,
                                                  "boundary velocity residual",
                                                  "A visualization output object that generates output for the velocity "
                                                  "residual at the top surface. The residual is computed at each point at the "
                                                  "surface as the difference between the modeled velocities and the input "
                                                  "data velocities for each vector component. The user has an option to choose "
                                                  "the input data as ascii data files (e.g. GPS velocities) with columns "
                                                  "in the same format as described for the 'ascii data' initial temperature plugin "
                                                  "or a velocity field computed from the GPlates program as described in the gplates "
                                                  "boundary velocity plugin."
                                                  "\n\n"
                                                  "Physical units: $\\frac{\\text{m}}{\\text{s}}$ or "
                                                  "$\\frac{\\text{m}}{\\text{year}}$, depending on settings in the input file.")
    }
  }
}
