/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>


namespace aspect
{
  namespace GeometryModel
  {
  template <int dim>
  Sphere<dim>::SphereGeometry::SphereGeometry()
  {}

  template <int dim>
  Point<dim>
  Sphere<dim>::SphereGeometry::
  push_forward(const Point<dim> &space_point) const
  {
	  /**
	   * The spherical_manifold push_forward function requires R_theta_co_phi,
	   * but the common interface is R_phi_theta. So we change it before
	   * we return it.
	   */
	  if(dim == 3)
	  {
		 const Point<dim> R_theta_co_phi (space_point(0),
				                           space_point(2),
										   space_point(1)+(90*(numbers::PI/180)));

	     return spherical_manifold.push_forward(R_theta_co_phi);
	  }
	  else
	  {
		  /**
		   * Following the deal.ii spherical_manifold convention,
		   * latitude is removed in 2d. TODO: It might be interesting
		   * to implement spacedim here.
		   */
			 const Point<dim> R_theta_co_phi (space_point(0),
					                           space_point(1));

		     return spherical_manifold.push_forward(R_theta_co_phi);
	  }

  }


  template <int dim>
  Point<dim>
  Sphere<dim>::SphereGeometry::
  pull_back(const Point<dim> &chart_point) const
  {
	  /**
	   * The spherical_manifold pull_back function returns R_theta_co_phi,
	   * but the common interface is R_phi_theta. So we change it before
	   * we return it.
	   * R, co_lat, long (0 to 2*pi) -> R, long (-pi to pi), lat
	   */
	  if(dim == 3)
	  {
			 const Point<dim> R_theta_co_phi = spherical_manifold.pull_back(chart_point);
		     return Point<dim>(R_theta_co_phi(0),
		    		           (R_theta_co_phi(2) > numbers::PI ? R_theta_co_phi(2) - 2 * numbers::PI :  R_theta_co_phi(2)),
							   (90*(numbers::PI/180))-R_theta_co_phi(1));
	  }
	  else
	  {
		  /**
		   * Following the deal.ii spherical_manifold convention,
		   * latitude is removed in 2d. TODO: It might be interesting
		   * to implement spacedim here.
		   */
			 const Point<dim> R_theta_co_phi = spherical_manifold.pull_back(chart_point);
		     return Point<dim>(R_theta_co_phi(0),
		    		           (R_theta_co_phi(1) > numbers::PI ? R_theta_co_phi(1) - 2 * numbers::PI :  R_theta_co_phi(1)));
	  }

  }


  template<int dim>
  const typename Interface<dim>::Manifold *
  Sphere<dim>::get_manifold() const
  {
    return &manifold;
  }
    template <int dim>
    void
    Sphere<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      GridGenerator::hyper_ball (coarse_grid,
                                 Point<dim>(),
                                 R);
      static const HyperBallBoundary<dim> boundary_ball(Point<dim>(), R);
      coarse_grid.set_boundary (0, boundary_ball);
    }


    template <int dim>
    std::set<types::boundary_id>
    Sphere<dim>::
    get_used_boundary_indicators () const
    {
      const types::boundary_id s[] = { 0 };
      return std::set<types::boundary_id>(&s[0],
                                          &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    Sphere<dim>::
    get_symbolic_boundary_names_map () const
    {
      static const std::pair<std::string,types::boundary_id> mapping("surface", 0);
      return std::map<std::string,types::boundary_id> (&mapping,
                                                       &mapping+1);
    }


    template <int dim>
    double
    Sphere<dim>::
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
    Sphere<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (R-position.norm(), 0.), maximal_depth());
    }



    template <int dim>
    Point<dim>
    Sphere<dim>::representative_point(const double depth) const
    {
      Point<dim> p;
      p(dim-1) = std::min (std::max(R - depth, 0.), maximal_depth());
      return p;
    }



    template <int dim>
    double
    Sphere<dim>::maximal_depth() const
    {
      return R;
    }

    template <int dim>
    double Sphere<dim>::radius () const
    {
      return R;
    }

    template <int dim>
    bool
    Sphere<dim>::has_curved_elements () const
    {
      return true;
    }



    template <int dim>
    bool
    Sphere<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(this->get_free_surface_boundary_indicators().size() == 0 ||
                  this->get_timestep_number() == 0,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(dynamic_cast<const InitialTopographyModel::ZeroTopography<dim>*>(&this->get_initial_topography_model()) != 0,
                  ExcMessage("After adding topography, this function can no longer be used to determine whether a point lies in the domain or not."));

      const double radius = point.norm();

      if (radius > R+std::numeric_limits<double>::epsilon()*R)
        return false;

      return true;
    }



    template <int dim>
    void
    Sphere<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Sphere");
        {
          prm.declare_entry ("Radius", "6371000",
                             Patterns::Double (0),
                             "Radius of the sphere. Units: m.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Sphere<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Sphere");
        {
          R            = prm.get_double ("Radius");
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
    ASPECT_REGISTER_GEOMETRY_MODEL(Sphere,
                                   "sphere",
                                   "Geometry model for sphere with a user specified radius. This geometry "
                                   "has only a single boundary, so the only valid boundary indicator to "
                                   "specify in the input file is ``0''. It can also be referenced by the "
                                   "symbolic name ``surface'' in input files.")
  }
}
