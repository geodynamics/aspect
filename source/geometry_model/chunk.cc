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


#include <aspect/geometry_model/chunk.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <iostream>
#include <fstream>



namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Chunk<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      std::vector<unsigned int> rep_vec(repetitions, repetitions+dim);
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 point1,
                                                 point2,
                                                 true);

      // Transform box into spherical chunk
      GridTools::transform (std_cxx11::bind(&SphericalManifold<dim>::push_forward,
                                            std_cxx11::cref(manifold),
                                            std_cxx11::_1),
                            coarse_grid);

      std::ofstream out ("grid-1.eps");
      GridOut grid_out;
      grid_out.write_eps (coarse_grid, out);

      // Deal with a curved mesh
      // Attach the real manifold to slot 15. we won't use it
      // during regular operation, but we set manifold_ids for all
      // cells, faces and edges immediately before refinement and
      // clear it again afterwards
      coarse_grid.set_manifold (15, manifold);

      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        cell->set_all_manifold_ids (15);

      // clear the manifold id from objects for which we have boundary
      // objects (and need boundary objects because at the time of
      // writing, only boundary objects provide normal vectors)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f))
            cell->face(f)->set_all_manifold_ids (numbers::invalid_manifold_id);

      for (typename Triangulation<dim>::active_cell_iterator
           cell = coarse_grid.begin_active();
           cell != coarse_grid.end(); ++cell)
            if (cell->at_boundary(0))
    		  std::cout << cell->face(0)->vertex(0) << std::endl;

      // deal.II wants boundary objects even for the straight boundaries
      // when using manifolds in the interior:
      static const StraightBoundary<dim> straight_boundary;
      coarse_grid.set_boundary	(4, straight_boundary);
      coarse_grid.set_boundary	(5, straight_boundary);

      // attach boundary objects to the curved boundaries:
      static const HyperShellBoundary<dim> boundary_shell;
      coarse_grid.set_boundary (0, boundary_shell);
      coarse_grid.set_boundary (1, boundary_shell);

      Point<dim> center;
      Point<dim> north, south;
      north[0] = 6371000.0 * std::cos(20.0*numbers::PI/180.0);
      south[0] = 6371000.0 * std::cos(45.0*numbers::PI/180.0);
      const double north_radius = std::sqrt(6371000.0*6371000.0+north[0]*north[0]);
      const double south_radius = std::sqrt(6371000.0*6371000.0+south[0]*south[0]);

      static const ConeBoundary<dim> boundary_cone_north(0,north_radius,center,north);
      coarse_grid.set_boundary (3, boundary_cone_north);
      static const ConeBoundary<dim> boundary_cone_south(0,south_radius,center,north);
            coarse_grid.set_boundary (2, boundary_cone_south);

      coarse_grid.signals.pre_refinement.connect (std_cxx11::bind (&set_manifold_ids,
                                                                   std_cxx11::ref(coarse_grid)));
      coarse_grid.signals.post_refinement.connect (std_cxx11::bind (&clear_manifold_ids,
                                                                    std_cxx11::ref(coarse_grid)));

    }



    template <int dim>
    std::set<types::boundary_id>
    Chunk<dim>::
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
    Chunk<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("inner",  0),
                  std::pair<std::string,types::boundary_id>("outer",  1),
                  std::pair<std::string,types::boundary_id>("west",   2),
                  std::pair<std::string,types::boundary_id>("east",   3)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }

          case 3:
          {
            static const std::pair<std::string,types::boundary_id> mapping[]
              = { std::pair<std::string,types::boundary_id>("inner",  0),
                  std::pair<std::string,types::boundary_id>("outer",  1),
                  std::pair<std::string,types::boundary_id>("west",   2),
                  std::pair<std::string,types::boundary_id>("east",   3),
                  std::pair<std::string,types::boundary_id>("south",  4),
                  std::pair<std::string,types::boundary_id>("north",  5)
                };

            return std::map<std::string,types::boundary_id> (&mapping[0],
                                                             &mapping[sizeof(mapping)/sizeof(mapping[0])]);
          }
        }

      Assert (false, ExcNotImplemented());
      return std::map<std::string,types::boundary_id>();
    }



    template <int dim>
    double
    Chunk<dim>::
    length_scale () const
    {
      // As described in the first ASPECT paper, a length scale of
      // 10km = 1e4m works well for the pressure scaling for earth
      // sized spherical shells. use a length scale that
      // yields this value for the R0,R1 corresponding to earth
      // but otherwise scales like (R1-R0)
      return 1e4 * maximal_depth() / (6336000.-3481000.);
    }



    template <int dim>
    double
    Chunk<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (point2[0]-position.norm(), 0.), maximal_depth());
    }



    template <int dim>
    Point<dim>
    Chunk<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // Choose a point at the mean longitude (and latitude)
      Point<dim> p = 0.5*(point2+point1);
      // at a depth beneath the top surface
      p[0] = point2[0]-depth;

      // Now convert to Cartesian coordinates
      return manifold.push_forward(p);
    }



    template <int dim>
    double
    Chunk<dim>::west_longitude () const
    {
    	if (dim == 2)
      return point1[1];
    	else
    		return point1[2];
    }



    template <int dim>
    double
    Chunk<dim>::east_longitude () const
    {
    	if (dim == 2)
      return point2[1];
    	else
    		return point2[2];
    }



    template <int dim>
    double
    Chunk<dim>::longitude_range () const
    {
    	if (dim == 2)
      return point2[1] - point1[1];
    	else
    		return point2[2] - point1[2];

    }



    template <int dim>
    double
    Chunk<dim>::south_latitude () const
    {
      if (dim == 3)
        return point1[1];
      else
        return 0;
    }



    template <int dim>
    double
    Chunk<dim>::north_latitude () const
    {
      if (dim==3)
        return point2[1];
      else
        return 0;
    }



    template <int dim>
    double
    Chunk<dim>::latitude_range () const
    {
      if (dim==3)
        return point2[1] - point1[1];
      else
        return 0;
    }



    template <int dim>
    double
    Chunk<dim>::maximal_depth() const
    {
      return point2[0]-point1[0];
    }



    template <int dim>
    double
    Chunk<dim>::inner_radius () const
    {
      return point1[0];
    }



    template <int dim>
    double
    Chunk<dim>::outer_radius () const
    {
      return point2[0];
    }



    template <int dim>
    bool
    Chunk<dim>::has_curved_elements() const
    {
      return true;
    }



    template <int dim>
    void
    Chunk<dim>::set_manifold_ids (Triangulation<dim> &triangulation)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end(); ++cell)
        cell->set_all_manifold_ids (15);
    }



    template <int dim>
    void
    Chunk<dim>::clear_manifold_ids (Triangulation<dim> &triangulation)
    {
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
            	cell->face(f)->set_all_manifold_ids (numbers::invalid_manifold_id);
        //cell->set_all_manifold_ids (numbers::invalid_manifold_id);
    }

    template <int dim>
    void
    Chunk<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {
          prm.declare_entry ("Chunk inner radius", "0",
                             Patterns::Double (0),
                             "Radius at the bottom surface of the chunk. Units: m.");
          prm.declare_entry ("Chunk outer radius", "1",
                             Patterns::Double (0),
                             "Radius at the top surface of the chunk. Units: m.");

          prm.declare_entry ("Chunk minimum longitude", "0",
                             Patterns::Double (-180, 360), // enables crossing of either hemisphere
                             "Minimum longitude of the chunk. Units: degrees.");
          prm.declare_entry ("Chunk maximum longitude", "1",
                             Patterns::Double (-180, 360), // enables crossing of either hemisphere
                             "Maximum longitude of the chunk. Units: degrees.");

          prm.declare_entry ("Chunk minimum latitude", "0",
                             Patterns::Double (-89, 89), // -90/90 degrees will crash because north/south boundary disappears
                             "Minimum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");
          prm.declare_entry ("Chunk maximum latitude", "1",
                             Patterns::Double (-89, 89), // -90/90 degrees will crash because north/south boundary disappears
                             "Maximum latitude of the chunk. This value is ignored "
                             "if the simulation is in 2d. Units: degrees.");

          prm.declare_entry ("Radius repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in radius.");
          prm.declare_entry ("Longitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in longitude.");
          prm.declare_entry ("Latitude repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in latitude. This value is ignored "
                             "if the simulation is in 2d");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Chunk<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {

          const double degtorad = dealii::numbers::PI/180;

          Assert (dim >= 2, ExcInternalError());
          Assert (dim <= 3, ExcInternalError());

          point1[0] = prm.get_double ("Chunk inner radius");
          point2[0] = prm.get_double ("Chunk outer radius");
          repetitions[0] = prm.get_integer ("Radius repetitions");

          if (dim == 2)
            {
              point1[1] = prm.get_double ("Chunk minimum longitude") * degtorad;
              point2[1] = prm.get_double ("Chunk maximum longitude") * degtorad;
              repetitions[1] = prm.get_integer ("Longitude repetitions");

              AssertThrow (point1[0] < point2[0],
                           ExcMessage ("Inner radius must be less than outer radius."));
              AssertThrow (point1[1] < point2[1],
                           ExcMessage ("Minimum longitude must be less than maximum longitude."));
              AssertThrow (point2[1] - point1[1] < 2.*numbers::PI,
                           ExcMessage ("Maximum - minimum longitude should be less than 360 degrees."));
            }

          // The pull back functions require radius,latitude,longitude order in 3D
          if (dim == 3)
            {
              point1[2] = prm.get_double ("Chunk minimum longitude") * degtorad;
              point2[2] = prm.get_double ("Chunk maximum longitude") * degtorad;
              repetitions[2] = prm.get_integer ("Longitude repetitions");
              point1[1] = 0.5*numbers::PI - prm.get_double ("Chunk minimum latitude") * degtorad;
              point2[1] = 0.5*numbers::PI - prm.get_double ("Chunk maximum latitude") * degtorad;
              repetitions[1] = prm.get_integer ("Latitude repetitions");

//              AssertThrow (point1[2] < point2[2],
//                           ExcMessage ("Minimum latitude must be less than maximum latitude."));
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
    ASPECT_REGISTER_GEOMETRY_MODEL(Chunk,
                                   "chunk",
                                   "A geometry which can be described as a chunk of a spherical shell, "
                                   "bounded by lines of longitude, latitude and radius. "
                                   "The minimum and maximum longitude, (latitude) and depth of the chunk "
                                   "is set in the parameter file. The chunk geometry labels its "
                                   "2*dim sides as follows: ``west'' and ``east'': minimum and maximum "
                                   "longitude, ``south'' and ``north'': minimum and maximum latitude, "
                                   "``inner'' and ``outer'': minimum and maximum radii. "
                                   "Names in the parameter files are as follows: "
                                   "Chunk (minimum || maximum) (longitude || latitude): "
                                   "edges of geographical quadrangle (in degrees)"
                                   "Chunk (inner || outer) radius: Radii at bottom and top of chunk"
                                   "(Longitude || Latitude || Radius) repetitions: "
                                   "number of cells in each coordinate direction.")
  }
}

