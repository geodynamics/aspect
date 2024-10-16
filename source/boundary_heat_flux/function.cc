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


#include <aspect/boundary_heat_flux/function.h>
#include <aspect/global.h>
#include <aspect/geometry_model/interface.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace BoundaryHeatFlux
  {
    template <int dim>
    std::vector<Tensor<1,dim>>
    Function<dim>::
    heat_flux (const types::boundary_id /*boundary_indicator*/,
               const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
               const MaterialModel::MaterialModelOutputs<dim> &/*material_model_outputs*/,
               const std::vector<Tensor<1,dim>> &normal_vectors) const
    {
      const unsigned int n_evaluation_points = material_model_inputs.n_evaluation_points();
      std::vector<Tensor<1,dim>> heat_flux(normal_vectors);

      for (unsigned int i=0; i<n_evaluation_points; ++i)
        {
          const Point<dim> position = material_model_inputs.position[i];
          if (coordinate_system == Utilities::Coordinates::cartesian)
            {
              heat_flux[i] *= boundary_heat_flux_function.value(position);
            }
          else if (coordinate_system == Utilities::Coordinates::spherical)
            {
              const std::array<double,dim> spherical_coordinates =
                aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
              Point<dim> point;

              for (unsigned int d=0; d<dim; ++d)
                point[d] = spherical_coordinates[d];

              heat_flux[i] *= boundary_heat_flux_function.value(point);
            }
          else if (coordinate_system == Utilities::Coordinates::depth)
            {
              const double depth = this->get_geometry_model().depth(position);
              Point<dim> point;
              point(0) = depth;

              heat_flux[i] *= boundary_heat_flux_function.value(point);
            }
          else
            {
              AssertThrow(false, ExcNotImplemented());
            }
        }

      return heat_flux;
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_heat_flux_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_heat_flux_function.set_time (this->get_time());
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary heat flux model");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry ("Coordinate system", "cartesian",
                             Patterns::Selection ("cartesian|spherical|depth"),
                             "A selection that determines the assumed coordinate "
                             "system for the function variables. Allowed values "
                             "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                             "are interpreted as r,phi or r,phi,theta in 2d/3d "
                             "respectively with theta being the polar angle. `depth' "
                             "will create a function, in which only the first "
                             "parameter is non-zero, which is interpreted to "
                             "be the depth of the point.");

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary heat flux model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }
        try
          {
            boundary_heat_flux_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary heat flux model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryHeatFlux
  {
    ASPECT_REGISTER_BOUNDARY_HEAT_FLUX_MODEL(Function,
                                             "function",
                                             "Implementation of a model in which the boundary "
                                             "heat flux is given in terms of an explicit formula "
                                             "that is elaborated in the parameters in section "
                                             "``Boundary heat flux model|Function''. The format of these "
                                             "functions follows the syntax understood by the "
                                             "muparser library, see "
                                             "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`."
                                             "\n\n"
                                             "The formula you describe in the mentioned "
                                             "section is a scalar value for the heat flux that is assumed "
                                             "to be the flux normal to the boundary, and that has the unit "
                                             "W/(m$^2$) (in 3d) or W/m (in 2d). Negative fluxes are "
                                             "interpreted as the flow of heat into the domain, and positive "
                                             "fluxes are interpreted as heat flowing out of the domain."
                                             "\n\n"
                                             "The symbol $t$ indicating time that "
                                             "may appear in the formulas for the prescribed "
                                             "heat flux is interpreted as having units "
                                             "seconds unless the global parameter ``Use "
                                             "years in output instead of seconds'' has "
                                             "been set.")
  }
}
