//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    SphericalShell<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      if (phi == 360)
        {
          GridGenerator::hyper_shell (coarse_grid,
                                      Point<dim>(),
                                      R0,
                                      R1,
                                      (dim==3) ? 96 : 12,
                                      true);
        }
      else if (phi == 90)
        {
          GridGenerator::quarter_hyper_shell (coarse_grid,
                                              Point<dim>(),
                                              R0,
                                              R1,0,
                                              true);
        }
      else if (phi == 180)
        {
          GridGenerator::half_hyper_shell (coarse_grid,
                                           Point<dim>(),
                                           R0,
                                           R1,0,
                                           true);
        }
      else
        {
          Assert (false, ExcInternalError());
        }

      static const HyperShellBoundary<dim> boundary;
      coarse_grid.set_boundary (0, boundary);
      coarse_grid.set_boundary (1, boundary);
    }


    template <int dim>
    std::set<unsigned char>
    SphericalShell<dim>::
    get_temperature_dirichlet_boundary_indicators () const
    {
      // see the documentation of the various GridGenerator::*hyper_shell
      // functions for a description of which boundary indicators are
      // set and how they correlate to what's used below
      const unsigned char s[] = { 0, 1 };
      return std::set<unsigned char>(&s[0],
                                     &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    std::set<unsigned char>
    SphericalShell<dim>::
    get_zero_velocity_boundary_indicators () const
    {
      // follow what is described in the documentation of this class.
      // see the documentation of the various GridGenerator::*hyper_shell
      // functions for a description of which boundary indicators are
      // set and how they correlate to what's used below
      if (phi == 360)
        {
          const unsigned char s[] = { 1 };
          return std::set<unsigned char>(&s[0],
                                         &s[sizeof(s)/sizeof(s[0])]);
        }
      else
        return std::set<unsigned char>();
    }



    template <int dim>
    std::set<unsigned char>
    SphericalShell<dim>::
    get_tangential_velocity_boundary_indicators () const
    {
      // follow what is described in the documentation of this class.
      // see the documentation of the various GridGenerator::*hyper_shell
      // functions for a description of which boundary indicators are
      // set and how they correlate to what's used below
      if (phi == 360)
        {
          const unsigned char s[] = { 0 };
          return std::set<unsigned char>(&s[0],
                                         &s[sizeof(s)/sizeof(s[0])]);
        }
      else
        {
          const unsigned char s[] = { 0, 1, 2, 3 };
          return std::set<unsigned char>(&s[0],
                                         &s[sizeof(s)/sizeof(s[0])]);
        }
    }


    template <int dim>
    double
    SphericalShell<dim>::
    length_scale () const
    {
      return 1e4;
    }


    template <int dim>
    double SphericalShell<dim>::inner_radius () const
    {
      return R0;
    }



    template <int dim>
    double SphericalShell<dim>::outer_radius () const
    {
      return R1;
    }



    template <int dim>
    double SphericalShell<dim>::opening_angle () const
    {
      return phi;
    }


    template <int dim>
    void
    SphericalShell<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Spherical shell");
        {
          prm.declare_entry ("Inner radius", "3481000",  // 6371-2890 in km
                             Patterns::Double (0),
                             "Inner radius of the spherical shell. Units: m.");
          prm.declare_entry ("Outer radius", "6336000",  // 6371-35 in km
                             Patterns::Double (0),
                             "Outer radius of the spherical shell. Units: m.");
          prm.declare_entry ("Opening angle", "360",
                             Patterns::Double (0, 360),
                             "Opening angle in degrees of the section of the shell "
			     "that we want to build. Units: degrees.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    SphericalShell<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Spherical shell");
        {
          R0            = prm.get_double ("Inner radius");
          R1            = prm.get_double ("Outer radius");
          phi = prm.get_double ("Opening angle");
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
  namespace GeometryModel
  {
    template class SphericalShell<deal_II_dimension>;

    ASPECT_REGISTER_GEOMETRY_MODEL(SphericalShell,
                                   "spherical shell",
                                   "A geometry representing a spherical shell or a pice of it. "
                                   "Inner and outer radii are read from the parameter file "
                                   "in subsection 'Spherical shell'.")
  }
}
