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


#include <aspect/geometry_model/cylinder.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Cylinder<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      GridGenerator::cylinder (coarse_grid,
                                 R, H/2.0);
      static const CylinderBoundary<dim> boundary_cylinder(R);
      coarse_grid.set_boundary (0, boundary_cylinder);
    }


    template <int dim>
    std::set<types::boundary_id>
    Cylinder<dim>::
    get_used_boundary_indicators () const
    {
      const types::boundary_id s[] = { 0 , 1, 2};
      return std::set<types::boundary_id>(&s[0],
                                          &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    double
    Cylinder<dim>::
    length_scale () const
    {
      return 1e4 * maximal_depth() / (6336000.-3481000.);
    }



    template <int dim>
    double
    Cylinder<dim>::depth(const Point<dim> &position) const
    {
      return 0.0; 
    }



    template <int dim>
    Point<dim>
    Cylinder<dim>::representative_point(const double depth) const
    {
      Point<dim> p;
      return p;
    }



    template <int dim>
    double
    Cylinder<dim>::maximal_depth() const
    {
      return H;
    }

    template <int dim>
    double Cylinder<dim>::radius () const
    {
      return R;
    }

    template <int dim>
    double Cylinder<dim>::height () const
    {
      return H;
    }

    template <int dim>
    void
    Cylinder<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Cylinder");
        {
          prm.declare_entry ("Radius", "1.0",
                             Patterns::Double (0),
                             "Radius of the cylinder. Units: m.");
          prm.declare_entry ("Height", "1.0",
                             Patterns::Double (0),
                             "Radius of the cylinder. Units: m.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Cylinder<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Cylinder");
        {
          R            = prm.get_double ("Radius");
          H            = prm.get_double ("Height");
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
    ASPECT_REGISTER_GEOMETRY_MODEL(Cylinder,
                                   "cylinder",
                                   "")
  }
}
