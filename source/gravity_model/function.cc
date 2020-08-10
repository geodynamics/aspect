/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/gravity_model/function.h>

#include <iostream>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Function<dim>::Function ()
      :
      function (dim)
    {}

    template <int dim>
    Tensor<1,dim>
    Function<dim>::
    gravity_vector (const Point<dim> &position) const
    {
      Tensor<1,dim> gravity;
      for (unsigned int d=0; d<dim; ++d)
        gravity[d] = function.value(position,d);
      return gravity;
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        function.set_time (this->get_time() / year_in_seconds);
      else
        function.set_time (this->get_time());
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Function");
        {
          Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Function");
        try
          {
            function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Gravity Model.Function'\n"
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
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(Function,
                                  "function",
                                  "Gravity is given in terms of an explicit formula "
                                  "that is elaborated in the parameters in section "
                                  "``Gravity model|Function''. The format of these "
                                  "functions follows the syntax understood by the "
                                  "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
