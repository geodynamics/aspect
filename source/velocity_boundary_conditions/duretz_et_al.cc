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


#include <aspect/velocity_boundary_conditions/duretz_et_al.h>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    namespace DuretzEtAl
    {
      namespace
      {
        //TODO: This is a copy of the code in duretz_et_al.cc. Unify it in
        // some way!
        // based on http://geodynamics.org/hg/cs/AMR/Discontinuous_Stokes with permission
        void _Inclusion(double pos[2], double r_inclusion, double eta, double *vx, double *vy, double *p)
        {
          const double min_eta = 1.0;
          const double max_eta = eta;
          const double epsilon = 1; //strain rate
          const double A(min_eta*(max_eta-min_eta)/(max_eta+min_eta));
          std::complex<double> phi, psi, dphi;
          const double offset[2]= {1.0, 1.0};
          double r2_inclusion = r_inclusion * r_inclusion;

          double x = pos[0]-offset[0];
          double y = pos[1]-offset[1];
          double r2 = x*x+y*y;

          std::complex<double> z(x,y);
          if (r2<r2_inclusion)
            {
              //inside the inclusion
              phi=0;
              dphi=0;
              psi=-4*epsilon*(max_eta*min_eta/(min_eta+max_eta))*z;
            }
          else
            {
              //outside the inclusion
              phi=-2*epsilon*A*r2_inclusion/z;
              dphi=-phi/z;
              psi=-2*epsilon*(min_eta*z+A*r2_inclusion*r2_inclusion/(z*z*z));
            }
          double visc = (r2<r2_inclusion)? max_eta : 1.0;
          std::complex<double> v = (phi - z*conj(dphi) - conj(psi))/(2.0*visc);
          *vx=v.real();
          *vy=v.imag();
          *p=-2*epsilon*dphi.real();
        }
      }



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
        _Inclusion(pos,0.2,eta_B, &velocity[0], &velocity[1], &pressure);

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
                                                   "benchmark.");
    }
  }
}
