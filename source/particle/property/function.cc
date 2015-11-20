/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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
        function (1)
      {}

      template <int dim>
      void
      Function<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                      const Vector<double> &,
                                                      const std::vector<Tensor<1,dim> > &,
                                                      std::vector<double> &data) const
      {
        data.push_back(function.value(position));
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      Function<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("function",1));
        return property_information;
      }


      template <int dim>
      void
      Function<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Function");
            {
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
          prm.enter_subsection("Tracers");
          {
            prm.enter_subsection("Function");
            try
              {
                function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Postprocess.Tracers.Function'\n"
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
                                        "Implementation of a model in which the tracer "
                                        "property is set by evaluating an explicit function "
                                        "at the initial position of each particle. The "
                                        "function is defined in the parameters in section "
                                        "``Tracers|Function''. The format of these "
                                        "functions follows the syntax understood by the "
                                        "muparser library, see Section~\\ref{sec:muparser-format}."
                                        "\n\n")
    }
  }
}

