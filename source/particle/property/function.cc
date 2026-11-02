// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/particle/property/function.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      Function<dim>::Function()
        :
        n_components (0)
      {}

      template <int dim>
      void
      Function<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                      std::vector<double> &data) const
      {
        // convert the position into the selected coordinate system
        const Utilities::NaturalCoordinate<dim> point = this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);

        for (unsigned int i = 0; i < n_components; ++i)
          {
            const double function_value = function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()), i);
            data.push_back(function_value);
          }
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      Function<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("function",n_components));
        return property_information;
      }


      template <int dim>
      void
      Function<dim>::declare_parameters (ParameterHandler &prm)
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
          prm.declare_entry ("Number of components", "1",
                             Patterns::Integer (0),
                             "The number of function components where each component is described "
                             "by a function expression delimited by a ';'.");
          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      Function<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Function");
        {
          n_components = prm.get_integer ("Number of components");
          try
            {
              function = std::make_unique<Functions::ParsedFunction<dim>>(n_components);
              function->parse_parameters (prm);
            }
          catch (...)
            {
              std::cerr << "ERROR: FunctionParser failed to parse\n"
                        << "\t'Particles.Function'\n"
                        << "with expression\n"
                        << "\t'" << prm.get("Function expression") << "'";
              throw;
            }

          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }
        prm.leave_subsection();
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(Function,
                                        "function",
                                        "Implementation of a model in which the particle "
                                        "property is set by evaluating an explicit function "
                                        "at the initial position of each particle. The "
                                        "function is defined in the parameters in section "
                                        "``Particles|Function''. The format of these "
                                        "functions follows the syntax understood by the "
                                        "muparser library, see "
                                        "{ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
    }
  }
}
