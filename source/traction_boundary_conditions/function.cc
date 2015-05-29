/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/traction_boundary_conditions/function.h>
#include <aspect/global.h>

namespace aspect
{
  namespace TractionBoundaryConditions
  {
    template <int dim>
    Function<dim>::Function ()
      :
      boundary_traction_function (dim)
    {}



    template <int dim>
    Tensor<1,dim>
    Function<dim>::
    traction (const Point<dim> &p,
              const Tensor<1,dim> &) const
    {
      Tensor<1,dim> traction;
      for (unsigned int d=0; d<dim; ++d)
        traction[d] = boundary_traction_function.value(p,d);

      return traction;
    }


    template <int dim>
    void
    Function<dim>::update()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        boundary_traction_function.set_time (this->get_time() / year_in_seconds);
      else
        boundary_traction_function.set_time (this->get_time());
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary traction model");
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
      prm.enter_subsection("Boundary traction model");
      {
        prm.enter_subsection("Function");
        try
          {
            boundary_traction_function.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Boundary traction model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
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
  namespace TractionBoundaryConditions
  {
    ASPECT_REGISTER_TRACTION_BOUNDARY_CONDITIONS(Function,
                                                 "function",
                                                 "Implementation of a model in which the boundary "
                                                 "traction is given in terms of an explicit formula "
                                                 "that is elaborated in the parameters in section "
                                                 "``Boundary traction model|Function''. "
                                                 "\n\n"
                                                 "The formula you describe in the mentioned "
                                                 "section is a semicolon separated list of traction components "
                                                 "for each of the $d$ components of the traction vector. "
                                                 "These $d$ formulas are interpreted as having units "
                                                 "Pa.")
  }
}
