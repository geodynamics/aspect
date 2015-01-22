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


#include <aspect/initial_conditions/spherical_shell.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <fstream>
#include <iostream>
#include <cstring>



namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    double
    SphericalHexagonalPerturbation<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // spherical shell. verify that it is indeed
      AssertThrow (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
                   (&this->get_geometry_model())
                   != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a spherical shell."));

      // this initial condition only makes sense if a boundary temperature
      // is prescribed. verify that it is indeed
      AssertThrow (&this->get_boundary_temperature()
                   != 0,
                   ExcMessage ("This initial condition can only be used if a boundary "
                               "temperature is prescribed."));

      const double R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                        (this->get_geometry_model()).outer_radius();

      // s = fraction of the way from
      // the inner to the outer
      // boundary; 0<=s<=1
      const double s = this->get_geometry_model().depth(position)
                       / this->get_geometry_model().maximal_depth();

      /* now compute an angular variation of the linear temperature field by
         stretching the variable s appropriately. note that the following
         formula leaves the end points s=0 and s=1 fixed, but stretches the
         region in between depending on the angle phi=atan2(x,y).

         For a plot, see
         http://www.wolframalpha.com/input/?i=plot+%28%282*sqrt%28x^2%2By^2%29-1%29%2B0.2*%282*sqrt%28x^2%2By^2%29-1%29*%281-%282*sqrt%28x^2%2By^2%29-1%29%29*sin%286*atan2%28x%2Cy%29%29%29%2C+x%3D-1+to+1%2C+y%3D-1+to+1
      */


      const double scale = ((dim==3)
                            ?
                            std::max(0.0,
                                     cos(3.14159 * fabs(position(2)/R1)))
                            :
                            1.0);
      const double phi   = std::atan2(position(0),position(1));
      const double s_mod = s
                           +
                           0.2 * s * (1-s) * std::sin(angular_mode*phi +(90 + 2*rotation_offset)*numbers::PI/180 ) * scale;

      return (this->get_boundary_temperature().maximal_temperature()*(s_mod)
              +
              this->get_boundary_temperature().minimal_temperature()*(1-s_mod));
    }



    template <int dim>
    void
    SphericalHexagonalPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Spherical hexagonal perturbation");
        {

          prm.declare_entry ("Angular mode", "6",
                             Patterns::Integer (),
                             "The number of convection cells to perturb the system with.");

          prm.declare_entry  ("Rotation offset", "-45",
                              Patterns::Double (),
                              "Amount of clockwise rotation in degrees to apply to "
                              "the perturbations. Default is set to -45 in order "
                              "to provide backwards compatibility.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    SphericalHexagonalPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Spherical hexagonal perturbation");
        {
          angular_mode = prm.get_integer ("Angular mode");
          rotation_offset = prm.get_double  ("Rotation offset");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    SphericalGaussianPerturbation<dim>::
    SphericalGaussianPerturbation()
    {

      /*       Note that the values we read in here have reasonable default values equation to
             the following:*/
      geotherm.resize(4);
      radial_position.resize(4);
      geotherm[0] = 1e0;
      geotherm[1] = 0.75057142857142856;
      geotherm[2] = 0.32199999999999995;
      geotherm[3] = 0.0;
      radial_position[0] =  0e0-1e-3;
      radial_position[1] =  0.16666666666666666;
      radial_position[2] =  0.83333333333333337;
      radial_position[3] =  1e0+1e-3;


    }

    template <int dim>
    double
    SphericalGaussianPerturbation<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // spherical shell. verify that it is indeed
      AssertThrow (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
                   (&this->get_geometry_model())
                   != 0,
                   ExcMessage ("This initial condition can only be used if the geometry "
                               "is a spherical shell."));

      // this initial condition only makes sense if a boundary temperature
      // is prescribed. verify that it is indeed
      AssertThrow (&this->get_boundary_temperature()
                   != 0,
                   ExcMessage ("This initial condition can only be used if a boundary "
                               "temperature is prescribed."));
      const double
      R0 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (this->get_geometry_model()).inner_radius(),
      R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (this->get_geometry_model()).outer_radius();
      const double dT = this->get_boundary_temperature().maximal_temperature()
                        - this->get_boundary_temperature().minimal_temperature();
      const double T0 = this->get_boundary_temperature().maximal_temperature()/dT;
      const double T1 = this->get_boundary_temperature().minimal_temperature()/dT;
      const double h = R1-R0;

      // s = fraction of the way from
      // the inner to the outer
      // boundary; 0<=s<=1
      const double r = position.norm();
      const double s = (r-R0)/h;

      const double scale=R1/(R1 - R0);
      const float eps = 1e-4;

      int indx = -1;
      for (unsigned int i=0; i<3; ++i)
        {
          if ((radial_position[i] - s) < eps && (radial_position[i+1] - s) > eps)
            {
              indx = i;
              break;
            }
        }
      Assert (indx >= 0, ExcInternalError());
      Assert (indx < 3,  ExcInternalError());
      int indx1 = indx + 1;
      const float dx = radial_position[indx1] - radial_position[indx];
      const float dy = geotherm[indx1] - geotherm[indx];

      const double InterpolVal = (( dx > 0.5*eps)
                                  ?
                                  // linear interpolation
                                  std::max(geotherm[3],geotherm[indx] + (s-radial_position[indx]) * (dy/dx))
                                  :
                                  // evaluate the point in the discontinuity
                                  0.5*( geotherm[indx] + geotherm[indx1] ));

      const double x = (scale - this->depth)*std::cos(angle);
      const double y = (scale - this->depth)*std::sin(angle);
      const double Perturbation = (sign * amplitude *
                                   std::exp( -( std::pow((position(0)*scale/R1-x),2)
                                                +
                                                std::pow((position(1)*scale/R1-y),2) ) / sigma));

      if (r > R1 - 1e-6*R1 || InterpolVal + Perturbation < T1)
        return T1*dT;
      else if (r < R0 + 1e-6*R0 || InterpolVal + Perturbation > T0 )
        return T0*dT;
      else
        return (InterpolVal + Perturbation)*dT;
    }



    template <int dim>
    void
    SphericalGaussianPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Spherical gaussian perturbation");
        {
          prm.declare_entry ("Angle", "0e0",
                             Patterns::Double (0),
                             "The angle where the center of the perturbation is placed.");
          prm.declare_entry ("Non-dimensional depth", "0.7",
                             Patterns::Double (0),
                             "The non-dimensional radial distance where the center of the "
                             "perturbation is placed.");
          prm.declare_entry ("Amplitude", "0.01",
                             Patterns::Double (0),
                             "The amplitude of the perturbation.");
          prm.declare_entry ("Sigma", "0.2",
                             Patterns::Double (0),
                             "The standard deviation of the Gaussian perturbation.");
          prm.declare_entry ("Sign", "1",
                             Patterns::Double (),
                             "The sign of the perturbation.");
          prm.declare_entry ("Filename for initial geotherm table", "initial-geotherm-table",
                             Patterns::FileName(),
                             "The file from which the initial geotherm table is to be read. "
                             "The format of the file is defined by what is read in "
                             "source/initial\\_conditions/spherical\\_shell.cc.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    SphericalGaussianPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Spherical gaussian perturbation");
        {
          angle = prm.get_double ("Angle");
          depth = prm.get_double ("Non-dimensional depth");
          amplitude  = prm.get_double ("Amplitude");
          sigma  = prm.get_double ("Sigma");
          sign  = prm.get_double ("Sign");
          initial_geotherm_table = prm.get ("Filename for initial geotherm table");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(SphericalHexagonalPerturbation,
                                       "spherical hexagonal perturbation",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed following an $N$-fold pattern in a specified "
                                       "direction from an otherwise spherically symmetric "
                                       "state. The class's name comes from previous versions "
                                       "when the only option was $N=6$.")

    ASPECT_REGISTER_INITIAL_CONDITIONS(SphericalGaussianPerturbation,
                                       "spherical gaussian perturbation",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed by a single Gaussian added to an "
                                       "otherwise spherically symmetric state. Additional "
                                       "parameters are read from the parameter file in subsection "
                                       "'Spherical gaussian perturbation'.")
  }
}
