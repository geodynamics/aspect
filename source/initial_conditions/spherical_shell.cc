//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/initial_conditions/spherical_shell.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/tensor.h>

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
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
              (this->geometry_model)
              != 0,
              ExcMessage ("This initial condition can only be used if the geometry "
                          "is a spherical shell."));

      const double
      R0 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (*this->geometry_model).inner_radius(),
      R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (*this->geometry_model).outer_radius();
      const double h = R1-R0;

      // s = fraction of the way from
      // the inner to the outer
      // boundary; 0<=s<=1
      const double r = position.norm();
      const double s = (r-R0)/h;

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
                                     cos(3.14159 * abs(position(2)/R1)))
                            :
                            1.0);
      const double phi   = std::atan2(position(0),position(1));
      const double s_mod = s
                           +
                           0.2 * s * (1-s) * std::sin(6*phi) * scale;

      return (this->boundary_temperature->maximal_temperature()*(1.0-s_mod)
              +
              this->boundary_temperature->minimal_temperature()*s_mod);
    }



    template <int dim>
    double
    SphericalGaussianPerturbation<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // this initial condition only makes sense if the geometry is a
      // spherical shell. verify that it is indeed
      Assert (dynamic_cast<const GeometryModel::SphericalShell<dim>*>
              (this->geometry_model)
              != 0,
              ExcMessage ("This initial condition can only be used if the geometry "
                          "is a spherical shell."));

      const double
      R0 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (*this->geometry_model).inner_radius(),
      R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&> (*this->geometry_model).outer_radius();
      const double h = R1-R0;

      // s = fraction of the way from
      // the inner to the outer
      // boundary; 0<=s<=1
      const double r = position.norm();
      const double s = (r-R0)/h;

      const double scale=R1/(R1 - R0);
      const float eps = 1e-4;

      const double geotherm[4] = { this->boundary_temperature->maximal_temperature(),
                                   this->boundary_temperature->maximal_temperature() - 1300,
                                   this->boundary_temperature->minimal_temperature() + 1200,
                                   this->boundary_temperature->minimal_temperature()
                                 };

      const double depth[4] = { R0-1e-2*R0,
                                R0+500e3,
                                R1-500e3,
                                R1+1e-2*R1
                              };

      int indx = -1;
      for (unsigned int i=0; i<3; ++i)
        {
          if ((depth[i] - r) < eps && (depth[i+1] - r) > eps)
            {
              indx = i;
              break;
            }
        }
      Assert (indx >= 0, ExcInternalError());
      Assert (indx < 3,  ExcInternalError());
      int indx1 = indx + 1;
      const float dx = depth[indx1] - depth[indx];
      const float dy = geotherm[indx1] - geotherm[indx];

      const double InterpolVal = (( dx > 0.5*eps)
                                  ?
                                  // linear interpolation
                                  std::max(geotherm[3],geotherm[indx] + (r-depth[indx]) * (dy/dx))
                                  :
                                  // evaluate the point in the discontinuity
                                  0.5*( geotherm[indx] + geotherm[indx1] ));

      const double x = (scale - this->depth)*std::cos(angle);
      const double y = (scale - this->depth)*std::sin(angle);
      const double Perturbation = (sign * amplitude * this->boundary_temperature->maximal_temperature() *
                                   std::exp( -( std::pow((position(0)*scale/R1-x),2)
                                                +
                                                std::pow((position(1)*scale/R1-y),2) ) / sigma));

      if (r > R1 - 1e-2*R1)
        return geotherm[3];
      else if (r < R0 + 1e-2*R0)
        return geotherm[0];
      else
        return InterpolVal + Perturbation;
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
                             "The radial distance where the center of the perturbation is placed.");
          prm.declare_entry ("Amplitude", "0.01",
                             Patterns::Double (0),
                             "The amplitude of the perturbation.");
          prm.declare_entry ("Sigma", "0.2",
                             Patterns::Double (0),
                             "The standard deviation of the Gaussian perturbation.");
          prm.declare_entry ("Sign", "1",
                             Patterns::Double (),
                             "The sign of the perturbation.");
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
    template class SphericalHexagonalPerturbation<deal_II_dimension>;
    ASPECT_REGISTER_INITIAL_CONDITIONS("spherical hexagonal perturbation", SphericalHexagonalPerturbation);

    template class SphericalGaussianPerturbation<deal_II_dimension>;
    ASPECT_REGISTER_INITIAL_CONDITIONS("spherical gaussian perturbation", SphericalGaussianPerturbation);
  }
}
