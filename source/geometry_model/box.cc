/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/simulator_signals.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Box<dim>::initialize ()
    {
      // Get pointer to initial topography model
      topo_model = const_cast<InitialTopographyModel::Interface<dim>*>(&this->get_initial_topography_model());
      // Check that initial topography is required.
      // If so, connect the initial topography function
      // to the right signals: It should be applied after
      // the final initial adaptive refinement and after a restart.
      if (Plugins::plugin_type_matches<InitialTopographyModel::ZeroTopography<dim>>(*topo_model) == false)
        {
          this->get_signals().pre_set_initial_state.connect(
            [&](typename parallel::distributed::Triangulation<dim> &tria)
          {
            this->topography(tria);
          }
          );
          this->get_signals().post_resume_load_user_data.connect(
            [&](typename parallel::distributed::Triangulation<dim> &tria)
          {
            this->topography(tria);
          }
          );
        }
    }


    template <int dim>
    void
    Box<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      std::vector<unsigned int> rep_vec(repetitions, repetitions+dim);
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 box_origin,
                                                 box_origin+extents,
                                                 true);

      // Tell p4est about the periodicity of the mesh.
      std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
      periodicity_vector;
      for (int i=0; i<dim; ++i)
        if (periodic[i])
          GridTools::collect_periodic_faces
          ( coarse_grid, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
            /*direction*/ i, periodicity_vector);

      if (periodicity_vector.size() > 0)
        coarse_grid.add_periodicity (periodicity_vector);
    }

    template <int dim>
    void
    Box<dim>::
    topography (typename parallel::distributed::Triangulation<dim> &grid) const
    {
      // Here we provide GridTools with the function to displace vertices
      // in the vertical direction by an amount specified by the initial topography model
      GridTools::transform(
        [&](const Point<dim> &p) -> Point<dim>
      {
        return this->add_topography(p);
      },
      grid);

      this->get_pcout() << "   Added initial topography to grid" << std::endl << std::endl;
    }


    template <int dim>
    Point<dim>
    Box<dim>::
    add_topography (const Point<dim> &x_y_z) const
    {
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; d++)
        surface_point[d] = x_y_z[d];

      // Get the surface topography at this point
      const double topo = topo_model->value(surface_point);

      // Compute the displacement of the z coordinate
      const double ztopo = (x_y_z[dim-1] - box_origin[dim-1]) / extents[dim-1] * topo;

      // Compute the new point
      Point<dim> x_y_ztopo = x_y_z;
      x_y_ztopo[dim-1] += ztopo;

      return x_y_ztopo;
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
    std::map<std::string,types::boundary_id>
    Box<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("left",   0),
                  std::pair<std::string,types::boundary_id>("right",  1),
                  std::pair<std::string,types::boundary_id>("bottom", 2),
                  std::pair<std::string,types::boundary_id>("top",    3)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("left",   0),
                  std::pair<std::string,types::boundary_id>("right",  1),
                  std::pair<std::string,types::boundary_id>("front",  2),
                  std::pair<std::string,types::boundary_id>("back",   3),
                  std::pair<std::string,types::boundary_id>("bottom", 4),
                  std::pair<std::string,types::boundary_id>("top",    5)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }


    template <int dim>
    std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
    Box<dim>::
    get_periodic_boundary_pairs () const
    {
      std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundaries;
      for ( unsigned int i=0; i<dim; ++i)
        if (periodic[i])
          periodic_boundaries.insert( std::make_pair( std::pair<types::boundary_id, types::boundary_id>(2*i, 2*i+1), i) );
      return periodic_boundaries;
    }

    template <int dim>
    Point<dim>
    Box<dim>::get_extents () const
    {
      return extents;
    }

    template <int dim>
    Point<dim>
    Box<dim>::get_origin () const
    {
      return box_origin;
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
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = position[d];

      // Get the surface topography at this point
      const double topo = topo_model->value(surface_point);

      const double d = extents[dim-1] + topo - (position(dim-1)-box_origin[dim-1]);
      return std::min (std::max (d, 0.), maximal_depth());
    }


    template <int dim>
    double
    Box<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return (position(dim-1)-box_origin[dim-1]) - extents[dim-1];
    }


    template <int dim>
    Point<dim>
    Box<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // choose a point on the center axis of the domain (without topography)
      Point<dim> p = extents/2+box_origin;
      p[dim-1] = extents[dim-1]+box_origin[dim-1]-depth;

      return p;
    }


    template <int dim>
    double
    Box<dim>::maximal_depth() const
    {
      return extents[dim-1] + topo_model->max_topography();
    }

    template <int dim>
    bool
    Box<dim>::has_curved_elements() const
    {
      return false;
    }

    template <int dim>
    bool
    Box<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(!this->get_parameters().mesh_deformation_enabled ||
                  this->simulator_is_past_initialization() == false,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()),
                  ExcMessage("After adding topography, this function can no longer be used "
                             "to determine whether a point lies in the domain or not."));

      for (unsigned int d = 0; d < dim; d++)
        if (point[d] > extents[d]+box_origin[d]+std::numeric_limits<double>::epsilon()*extents[d] ||
            point[d] < box_origin[d]-std::numeric_limits<double>::epsilon()*extents[d])
          return false;

      return true;
    }

    template <int dim>
    std::array<double,dim>
    Box<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; i++)
        position_array[i] = position_point(i);

      return position_array;
    }


    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    Box<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }


    template <int dim>
    Point<dim>
    Box<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; i++)
        position_point[i] = position_tensor[i];

      return position_point;
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
          prm.declare_entry ("X extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in x-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Y extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in y-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Z extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in z-direction. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("Box origin X coordinate", "0.",
                             Patterns::Double (),
                             "X coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Y coordinate", "0.",
                             Patterns::Double (),
                             "Y coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Z coordinate", "0.",
                             Patterns::Double (),
                             "Z coordinate of box origin. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction.");

          prm.declare_entry ("X periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in X direction");
          prm.declare_entry ("Y periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Y direction");
          prm.declare_entry ("Z periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Z direction");

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
          box_origin[0] = prm.get_double ("Box origin X coordinate");
          extents[0] = prm.get_double ("X extent");
          periodic[0] = prm.get_bool ("X periodic");
          repetitions[0] = prm.get_integer ("X repetitions");

          if (dim >= 2)
            {
              box_origin[1] = prm.get_double ("Box origin Y coordinate");
              extents[1] = prm.get_double ("Y extent");
              periodic[1] = prm.get_bool ("Y periodic");
              repetitions[1] = prm.get_integer ("Y repetitions");
            }

          if (dim >= 3)
            {
              // Use dim-1 instead of 2 to avoid compiler warning in 2d:
              box_origin[dim-1] = prm.get_double ("Box origin Z coordinate");
              extents[dim-1] = prm.get_double ("Z extent");
              periodic[dim-1] = prm.get_bool ("Z periodic");
              repetitions[dim-1] = prm.get_integer ("Z repetitions");
            }
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
                                   "and top boundaries (see also the documentation of the deal.II class "
                                   "``GeometryInfo''). You can also use symbolic names ``left'', ``right'', "
                                   "etc., to refer to these boundaries in input files. "
                                   "It is also possible to add initial topography to the box model. Note however that "
                                   "this is done after the last initial adaptive refinement cycle. "
                                   "Also, initial topography is supposed to be small, as it is not taken into account "
                                   "when depth or a representative point is computed. ")


  }
}
