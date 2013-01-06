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


#include <aspect/geometry_model/box.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Box<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      GridGenerator::hyper_rectangle (coarse_grid,
                                      Point<dim>(),
                                      extents);
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        coarse_grid.begin_active()->face(f)->set_boundary_indicator(f);
    }


    template <int dim>
    std::set<types::boundary_id>
    Box<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }


    template <int dim>
    Point<dim>
    Box<dim>::get_extents () const
    {
      return extents;
    }


    template <int dim>
    double
    Box<dim>::
    length_scale () const
    {
      return 0.01*extents[0];
    }


    template <int dim>
    double
    Box<dim>::depth(const Point<dim> &position) const
    {
      const double d = maximal_depth()-position(dim-1);

      Assert (d >= 0, ExcInternalError());
      Assert (d <= maximal_depth(), ExcInternalError());

      return d;
    }


    template <int dim>
    Point<dim>
    Box<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // choose a point on the center axis of the domain
      Point<dim> p = extents/2;
      p[dim-1] = maximal_depth() - depth;
      return p;
    }


    template <int dim>
    double
    Box<dim>::maximal_depth() const
    {
      return extents[dim-1];
    }


    template <int dim>
    void
    Box<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("X extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in x-direction. Units: m.");
          prm.declare_entry ("Y extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in y-direction. Units: m.");
          prm.declare_entry ("Z extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in z-direction. This value is ignored "
                             "if the simulation is in 2d Units: m.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Box<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          extents[0] = prm.get_double ("X extent");

          if (dim >= 2)
            extents[1] = prm.get_double ("Y extent");

          if (dim >= 3)
            extents[2] = prm.get_double ("Z extent");
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
    ASPECT_REGISTER_GEOMETRY_MODEL(Box,
                                   "box",
                                   "A box geometry parallel to the coordinate directions. "
                                   "The extent of the box in each coordinate direction "
                                   "is set in the parameter file. The box geometry labels its "
                                   "2*dim sides as follows: in 2d, boundary indicators 0 through 3 "
                                   "denote the left, right, bottom and top boundaries; in 3d, boundary "
                                   "indicators 0 through 5 indicate left, right, front, back, bottom "
                                   "and top boundaries. See also the documentation of the deal.II class "
                                   "``GeometryInfo''.")
  }
}
