//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/initial_conditions_spherical_shell.h>

#include <deal.II/base/tensor.h>

namespace aspect
{
  namespace InitialConditions
  {
    template <int dim>
    double
    SphericalShellPerturbed<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      return 0;
    }

    template <int dim>
    void
    SphericalShellPerturbed<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Spherical shell perturbed");
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
      prm.declare_entry ("Gaussian perturbation", "true",
                         Patterns::Bool (),
                         "The sign of the perturbation.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    SphericalShellPerturbed<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Spherical shell perturbed");
        {
	  angle = prm.get_double ("Angle");
	  depth = prm.get_double ("Non-dimensional depth");
	  amplitude  = prm.get_double ("Amplitude");
	  sigma  = prm.get_double ("Sigma");
	  sign  = prm.get_double ("Sign");
	  gaussian_perturbation  = prm.get_bool ("Gaussian perturbation");
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
    template class SphericalShellPerturbed<deal_II_dimension>;
    ASPECT_REGISTER_INITIAL_CONDITIONS("Spherical shell perturbed", SphericalShellPerturbed);
  }
}
