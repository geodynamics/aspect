/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <aspect/particle/property/function.h>

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
        for (unsigned int i = 0; i < n_components; i++)
          data.push_back(function->value(position, i));
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      Function<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("function",n_components));
        return property_information;
      }


      template <int dim>
      void
      Function<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Function");
            {
              prm.declare_entry ("Number of components", "1",
                                 Patterns::Integer (0),
                                 "The number of function components where each component is described "
                                 "by a function expression delimited by a ';'.");
              Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      Function<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Function");
            n_components = prm.get_integer ("Number of components");
            try
              {
                function.reset (new Functions::ParsedFunction<dim>(n_components));
                function->parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Postprocess.Particles.Function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'";
                throw;
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
                                        "muparser library, see Section~\\ref{sec:muparser-format}.")
    }
  }
}

