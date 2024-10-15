/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/postprocess/visualization/surface_dynamic_topography.h>
#include <aspect/postprocess/dynamic_topography.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      SurfaceDynamicTopography<dim>::
      SurfaceDynamicTopography ()
        :
        DataPostprocessorScalar<dim> ("surface_dynamic_topography",
                                      update_quadrature_points),
        Interface<dim>("m")
      {}



      template <int dim>
      void
      SurfaceDynamicTopography<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        // Initialize everything to zero, so that we can ignore faces we are
        // not interested in (namely, those not labeled as 'top' or 'bottom')
        for (auto &quantity : computed_quantities)
          quantity(0) = 0;

        const Postprocess::DynamicTopography<dim> &dynamic_topography =
          this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

        auto cell = input_data.template get_cell<dim>();

        // We only want to output dynamic topography at the top and bottom
        // boundary, so only compute it if the current cell has
        // a face at the top or bottom boundary.
        bool cell_at_top_or_bottom_boundary = false;
        for (const unsigned int f : cell->face_indices())
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top" ||
               this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "bottom"))
            cell_at_top_or_bottom_boundary = true;

        if (cell_at_top_or_bottom_boundary)
          {
            std::vector<Point<dim>> quadrature_points(input_data.evaluation_points.size());
            for (unsigned int i=0; i<input_data.evaluation_points.size(); ++i)
              quadrature_points[i] = this->get_mapping().transform_real_to_unit_cell(cell,input_data.evaluation_points[i]);

            const Quadrature<dim> quadrature_formula(quadrature_points);

            FEValues<dim> fe_volume_values (this->get_mapping(),
                                            this->get_fe(),
                                            quadrature_formula,
                                            update_values);

            fe_volume_values.reinit(cell);

            std::vector<double> dynamic_topography_values(quadrature_formula.size());

            // It might seem unintuitive to use the extractor for the temperature block,
            // but that is where dynamic_topography.topography_vector() stores the values.
            // See the documentation of that function for more details.
            fe_volume_values[this->introspection().extractors.temperature].get_function_values(dynamic_topography.topography_vector(),
                dynamic_topography_values);

            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              computed_quantities[q](0) = dynamic_topography_values[q];
          }
      }



      /**
       * Register the other postprocessor that we need: DynamicTopography
       */
      template <int dim>
      std::list<std::string>
      SurfaceDynamicTopography<dim>::required_other_postprocessors() const
      {
        return {"dynamic topography"};
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceDynamicTopography,
                                                  "surface dynamic topography",
                                                  "A visualization output object that generates output "
                                                  "for the dynamic topography at the top and bottom of the model space. "
                                                  "The actual computation of this topography is handled inside the "
                                                  "'dynamic topography' postprocessor, please check its documentation "
                                                  "for details about the numerical methods."
                                                  "\n\n"
                                                  "In contrast to the `dynamic topography' visualization postprocessor, this "
                                                  "plugin really only evaluates the dynamic topography at faces of cells "
                                                  "that are adjacent to `bottom' and `top' boundaries, and only outputs "
                                                  "information on the surface of the domain, rather than padding the "
                                                  "information with zeros in the interior of the domain."
                                                  "\n\n"
                                                  "Physical units: \\si{\\meter}.")
    }
  }
}
