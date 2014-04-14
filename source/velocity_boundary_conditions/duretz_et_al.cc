/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/velocity_boundary_conditions/duretz_et_al.h>


namespace aspect
{
  namespace DuretzEtAl
  {
    namespace AnalyticSolutions
    {
      // this function is defined in source/postprocessor/duretz_et_al.cc
      void _Inclusion(double pos[2], double r_inclusion, double eta, double *vx, double *vy, double *p);
    }
  }


  namespace VelocityBoundaryConditions
  {
    namespace DuretzEtAl
    {
      template <int dim>
      Inclusion<dim>::Inclusion ()
        :
        eta_B (1e3)
      {}



      template <int dim>
      Tensor<1,dim>
      Inclusion<dim>::
      boundary_velocity (const Point<dim> &p) const
      {
        Assert (dim == 2, ExcNotImplemented());

        double pos[2]= {p(0),p(1)};

        Tensor<1,dim> velocity;
        double pressure;
        aspect::DuretzEtAl::AnalyticSolutions::_Inclusion
        (pos,0.2,eta_B, &velocity[0], &velocity[1], &pressure);

        return velocity;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    namespace DuretzEtAl
    {
      ASPECT_REGISTER_VELOCITY_BOUNDARY_CONDITIONS(Inclusion,
                                                   "inclusion",
                                                   "Implementation of the velocity boundary conditions for the "
                                                   "``inclusion'' benchmark. See the manual and the Kronbichler, Heister "
                                                   "and Bangerth paper on ASPECT for more information about this "
                                                   "benchmark.")
    }
  }
}
