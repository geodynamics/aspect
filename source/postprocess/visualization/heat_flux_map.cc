/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/heat_flux_map.h>
#include <aspect/postprocess/heat_flux_map.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/boundary_velocity/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      HeatFluxMap<dim>::
      HeatFluxMap ()
        :
        DataPostprocessorScalar<dim> ("heat_flux_map",
                                      update_quadrature_points)
      {}



      template <int dim>
      void
      HeatFluxMap<dim>::update ()
      {
        if (output_point_wise_heat_flux)
          heat_flux_density_solution = Postprocess::internal::compute_dirichlet_boundary_heat_flux_solution_vector(*this);
        else
          heat_flux_and_area = internal::compute_heat_flux_through_boundary_faces (*this);
      }



      template <int dim>
      void
      HeatFluxMap<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        for (unsigned int q=0; q<computed_quantities.size(); ++q)
          computed_quantities[q](0) = 0;

#if DEAL_II_VERSION_GTE(9,3,0)
        auto cell = input_data.template get_cell<dim>();
#else
        auto cell = input_data.template get_cell<DoFHandler<dim> >();
#endif

        if (output_point_wise_heat_flux)
          {
            bool cell_at_top_or_bottom_boundary = false;
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
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

                std::vector<double> heat_flux_values(quadrature_formula.size());
                fe_volume_values[this->introspection().extractors.temperature].get_function_values(heat_flux_density_solution, heat_flux_values);

                for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                  computed_quantities[q](0) = heat_flux_values[q];
              }
          }
        else
          {
            double heat_flux = 0.0;

            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->at_boundary(f) &&
                  (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top" ||
                   this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "bottom"))
                {
                  // add heatflow for this face
                  heat_flux += heat_flux_and_area[cell->active_cell_index()][f].first /
                               heat_flux_and_area[cell->active_cell_index()][f].second;
                }

            for (unsigned int q=0; q<computed_quantities.size(); ++q)
              computed_quantities[q](0) = heat_flux;
          }
      }



      template <int dim>
      void
      HeatFluxMap<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Heat flux map");
            {
              prm.declare_entry("Output point wise heat flux",
                                "false",
                                Patterns::Bool(),
                                "A boolean flag that controls whether to output the "
                                "heat flux map as a point wise value, or as a cell-wise "
                                "averaged value. The point wise output is more "
                                "accurate, but it currently omits prescribed heat "
                                "flux values at boundaries and advective heat flux "
                                "that is caused by velocities non-tangential to boundaries. "
                                "If you do not use these two features it is recommended "
                                "to switch this setting on to benefit from the increased "
                                "output resolution.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      HeatFluxMap<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Heat flux map");
            {
              output_point_wise_heat_flux = prm.get_bool("Output point wise heat flux");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(HeatFluxMap,
                                                  "heat flux map",
                                                  "A visualization output object that generates output for "
                                                  "the heat flux density across the top and bottom boundary "
                                                  "in outward direction. "
                                                  "The heat flux is computed as sum "
                                                  "of advective heat flux and conductive heat "
                                                  "flux through Neumann boundaries, both "
                                                  "computed as integral over the boundary area, "
                                                  "and conductive heat flux through Dirichlet "
                                                  "boundaries, which is computed using the "
                                                  "consistent boundary flux method as described "
                                                  "in ``Gresho, Lee, Sani, Maslanik, Eaton (1987). "
                                                  "The consistent Galerkin FEM for computing "
                                                  "derived boundary quantities in thermal and or "
                                                  "fluids problems. International Journal for "
                                                  "Numerical Methods in Fluids, 7(4), 371-394.'' "
                                                  "If only conductive heat flux through Dirichlet "
                                                  "boundaries is of interest, the "
                                                  "postprocessor can produce output of higher resolution "
                                                  "by evaluating the CBF solution vector point-wise "
                                                  "instead of computing cell-wise averaged values.")
    }
  }
}
