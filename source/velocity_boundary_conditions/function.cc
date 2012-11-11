/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */


#include <aspect/velocity_boundary_conditions/function.h>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    template <int dim>
    Function<dim>::Function ()
      :
      boundary_velocity_function (dim)
    {}



    template <int dim>
    Tensor<1,dim>
    Function<dim>::
    boundary_velocity (const Point<dim> &p) const
    {
      Tensor<1,dim> velocity;
      for (unsigned int d=0; d<dim; ++d)
        velocity[d] = boundary_velocity_function.value(p,d);

      return velocity;
    }


    template <int dim>
    void
    Function<dim>::set_current_time (const double time)
    {
      boundary_velocity_function.set_time(time);
    }


    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary velocity model");
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
      prm.enter_subsection("Boundary velocity model");
      {
        prm.enter_subsection("Function");
        {
          boundary_velocity_function.parse_parameters (prm);
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
  namespace VelocityBoundaryConditions
  {
    ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(Function,
                                                 "function",
                                                 "Implementation of a model in which the boundary "
                                                 "velocity is given in terms of an explicit formula "
                                                 "that is elaborated in the parameters in section "
                                                 "``Boundary velocity model|Function''.")
  }
}
