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


#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    SphericalShell<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      AssertThrow (phi == 360 || phi == 90 || dim!=3, ExcNotImplemented());

      if (phi == 360)
        {
          GridGenerator::hyper_shell (coarse_grid,
                                      Point<dim>(),
                                      R0,
                                      R1,
                                      (n_cells_along_circumference == 0
                                       ?
                                       // automatic choice that leads to reasonable
                                       // meshes with the typical aspect ratio of
                                       // the Earth
                                       (dim==3 ? 96 : 12)
                                       :
                                       // user choice
                                       n_cells_along_circumference),
                                      true);
        }
      else if (phi == 90)
        {
          GridGenerator::quarter_hyper_shell (coarse_grid,
                                              Point<dim>(),
                                              R0,
                                              R1,
                                              0,
                                              true);
        }
      else if (phi == 180)
        {
          GridGenerator::half_hyper_shell (coarse_grid,
                                           Point<dim>(),
                                           R0,
                                           R1,
                                           0,
                                           true);
        }
      else
        {
          Assert (false, ExcInternalError());
        }

      // Use a manifold description for all cells. use manifold_id 99 in order
      // not to step on the boundary indicators used below
      static const SphericalManifold<dim> spherical_manifold;
      coarse_grid.set_manifold (99, spherical_manifold);

      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        cell->set_all_manifold_ids (99);

      // clear the manifold id from objects for which we have boundary
      // objects (and need boundary objects because at the time of
      // writing, only boundary objects provide normal vectors)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f))
            cell->face(f)->set_all_manifold_ids (numbers::invalid_manifold_id);

      // deal.II wants boundary objects even for the straight boundaries
      // when using manifolds in the interior:
      static const StraightBoundary<dim> straight_boundary;
      std::set<types::boundary_id> ids = get_used_boundary_indicators();
      for (std::set<types::boundary_id>::iterator it = ids.begin();
           it!=ids.end(); ++it)
        if (*it > 1)
          coarse_grid.set_boundary (*it, straight_boundary);

      // attach boundary objects to the curved boundaries:
      static const HyperShellBoundary<dim> boundary_shell;
      coarse_grid.set_boundary (0, boundary_shell);
      coarse_grid.set_boundary (1, boundary_shell);
    }


    template <int dim>
    std::set<types::boundary_id>
    SphericalShell<dim>::
    get_used_boundary_indicators () const
    {
      // follow what is described in the documentation of this class.
      // see the documentation of the various GridGenerator::*hyper_shell
      // functions for a description of which boundary indicators are
      // set and how they correlate to what's used below
      if (phi == 360)
        {
          const types::boundary_id s[] = { 0, 1 };
          return std::set<types::boundary_id>(&s[0],
                                              &s[sizeof(s)/sizeof(s[0])]);
        }
      else if (phi == 90 && dim == 3)
        {
          const types::boundary_id s[] = { 0, 1, 2, 3, 4};
          return std::set<types::boundary_id>(&s[0],
                                              &s[sizeof(s)/sizeof(s[0])]);
        }
      else
        {
          const types::boundary_id s[] = { 0, 1, 2, 3 };
          return std::set<types::boundary_id>(&s[0],
                                              &s[sizeof(s)/sizeof(s[0])]);
        }
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    SphericalShell<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id> ("inner", 0),
                  std::pair<std::string,types::boundary_id> ("outer", 1),
                  std::pair<std::string,types::boundary_id> ("left",  2),
                  std::pair<std::string,types::boundary_id> ("right", 3)
                };

            if (phi == 360)
              return std::map<std::string,types::boundary_id> (&mapping[0],
                                                               &mapping[2]);
            else
              return std::map<std::string,types::boundary_id> (&mapping[0],
                                                               &mapping[4]);
          }

          case 3:
          {
            if (phi == 360)
              {
                static const std::pair<std::string,types::boundary_id> mapping[]
                  = { std::pair<std::string,types::boundary_id>("inner", 0),
                      std::pair<std::string,types::boundary_id>("outer", 1)
                    };

                return std::map<std::string,types::boundary_id> (&mapping[0],
                                                                 &mapping[2]);
              }
            else if (phi == 90)
              {
                static const std::pair<std::string,types::boundary_id> mapping[]
                  = { std::pair<std::string,types::boundary_id>("inner", 0),
                      std::pair<std::string,types::boundary_id>("outer", 1),
                      std::pair<std::string,types::boundary_id>("east",  2),
                      std::pair<std::string,types::boundary_id>("west",  3),
                      std::pair<std::string,types::boundary_id>("south", 4)
                    };

                return std::map<std::string,types::boundary_id> (&mapping[0],
                                                                 &mapping[5]);
              }
            else
              Assert (false, ExcNotImplemented());
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    double
    SphericalShell<dim>::
    length_scale () const
    {
      // as described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. use a length scale that
      // yields this value for the R0,R1 corresponding to earth
      // but otherwise scales like (R1-R0)
      return 1e4 * maximal_depth() / (6336000.-3481000.);
    }



    template <int dim>
    double
    SphericalShell<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (R1-position.norm(), 0.), maximal_depth());
    }



    template <int dim>
    Point<dim>
    SphericalShell<dim>::representative_point(const double depth) const
    {
      Point<dim> p;
      p(dim-1) = std::min (std::max(R1 - depth, R0), R1);
      return p;
    }



    template <int dim>
    double
    SphericalShell<dim>::maximal_depth() const
    {
      return R1-R0;
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
    bool
    SphericalShell<dim>::has_curved_elements () const
    {
      return true;
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

          prm.declare_entry ("Cells along circumference", "0",
                             Patterns::Integer (0),
                             "The number of cells in circumferential direction that are "
                             "created in the coarse mesh in 2d. If zero, this number "
                             "is chosen automatically in a way that produces meshes "
                             "in which cells have a reasonable aspect ratio for models "
                             "in which the depth of the mantle is roughly that of the "
                             "Earth. For planets with much shallower mantles and larger "
                             "cores, you may want to chose a larger number to avoid "
                             "cells that are elongated in tangential and compressed in "
                             "radial direction."
                             "\n\n"
                             "In 3d, the number of cells is computed differently and does "
                             "not have an easy interpretation. Valid values for this parameter "
                             "in 3d are 0 (let this class choose), 6, 12 and 96. "
                             "Other possible values may be discussed in the documentation "
                             "of the deal.II function GridGenerator::hyper_shell. "
                             "The parameter is best left at its default in 3d."
                             "\n\n"
                             "In either case, this parameter is ignored unless the opening "
                             "angle of the domain is 360 degrees.");
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
          R0  = prm.get_double ("Inner radius");
          R1  = prm.get_double ("Outer radius");
          phi = prm.get_double ("Opening angle");
          n_cells_along_circumference = prm.get_integer ("Cells along circumference");
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
    ASPECT_REGISTER_GEOMETRY_MODEL(SphericalShell,
                                   "spherical shell",
                                   "A geometry representing a spherical shell or a pice of it. "
                                   "Inner and outer radii are read from the parameter file "
                                   "in subsection 'Spherical shell'."
                                   "\n\n"
                                   "The model assigns boundary indicators as follows: In 2d, "
                                   "inner and outer boundaries get boundary indicators zero "
                                   "and one, and if the opening angle set in the input file "
                                   "is less than 360, then left and right boundaries are "
                                   "assigned indicators two and three. These boundaries can "
                                   "also be referenced using the symbolic names 'inner', 'outer' "
                                   "and (if applicable) 'left', 'right'."
                                   "\n\n"
                                   "In 3d, inner and "
                                   "outer indicators are treated as in 2d. If the opening "
                                   "angle is chosen as 90 degrees, i.e., the domain is the "
                                   "intersection of a spherical shell and the first octant, "
                                   "then indicator 2 is at the face $x=0$, 3 at $y=0$, "
                                   "and 4 at $z=0$. These last three boundaries can then also "
                                   "be referred to as 'east', 'west' and 'south' symbolically "
                                   "in input files.")
  }
}
