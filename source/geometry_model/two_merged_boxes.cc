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


#include <aspect/geometry_model/two_merged_boxes.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/base/std_cxx1x/bind.h>

namespace aspect
{
  namespace GeometryModel
  {

    template <int dim>
    void
    TwoMergedBoxes<dim>::
    set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      // iterate over all active cells and (re)set the boundary indicators
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active();
           cell != triangulation.end();
           ++cell)
        {

          // first set the default boundary indicators
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->face(f)->at_boundary())
#if DEAL_II_VERSION_GTE(8,3,0)
              cell->face(f)->set_boundary_id (f);
#else
              cell->face(f)->set_boundary_indicator (f);
#endif

          if (cell->face(0)->at_boundary())
            // set the lithospheric part of the left boundary to indicator 2*dim
            if (cell->face(0)->vertex(GeometryInfo<dim-1>::vertices_per_cell-1)[dim-1] > height_lith)
#if DEAL_II_VERSION_GTE(8,3,0)
              cell->face(0)->set_boundary_id (2*dim);
#else
              cell->face(0)->set_boundary_indicator (2*dim);
#endif

          if (cell->face(1)->at_boundary())
            // set the lithospheric part of the right boundary to indicator 2*dim+1
            if (cell->face(1)->vertex(GeometryInfo<dim-1>::vertices_per_cell-1)[dim-1] > height_lith)
#if DEAL_II_VERSION_GTE(8,3,0)
              cell->face(1)->set_boundary_id (2*dim+1);
#else
              cell->face(1)->set_boundary_indicator (2*dim+1);
#endif

          if (dim==3)
            {
              // set the lithospheric part of the front boundary to indicator 2*dim+2
              if (cell->face(2)->at_boundary())
                if (cell->face(2)->vertex(GeometryInfo<dim-1>::vertices_per_cell-1)[dim-1] > height_lith)
#if DEAL_II_VERSION_GTE(8,3,0)
                  cell->face(2)->set_boundary_id (2*dim+2);
#else
                  cell->face(2)->set_boundary_indicator (2*dim+2);
#endif

              // set the lithospheric part of the back boundary to indicator 2*dim+3
              if (cell->face(3)->at_boundary())
                if (cell->face(3)->vertex(GeometryInfo<dim-1>::vertices_per_cell-1)[dim-1] > height_lith)
#if DEAL_II_VERSION_GTE(8,3,0)
                  cell->face(3)->set_boundary_id (2*dim+3);
#else
                  cell->face(3)->set_boundary_indicator (2*dim+3);
#endif
            }
        }
    }

    template <int dim>
    void
    TwoMergedBoxes<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &total_coarse_grid) const
    {
      std::vector<unsigned int> lower_rep_vec(lower_repetitions, lower_repetitions+dim);
      std::vector<unsigned int> upper_rep_vec(upper_repetitions, upper_repetitions+dim);

      // the two triangulations that will be merged
      Triangulation<dim> lower_coarse_grid;
      Triangulation<dim> upper_coarse_grid;

      // create lower_coarse_grid mesh
      GridGenerator::subdivided_hyper_rectangle (lower_coarse_grid,
                                                 lower_rep_vec,
                                                 lower_box_origin,
                                                 lower_box_origin+lower_extents,
                                                 false);

      // create upper_coarse_grid mesh
      GridGenerator::subdivided_hyper_rectangle (upper_coarse_grid,
                                                 upper_rep_vec,
                                                 upper_box_origin,
                                                 upper_box_origin+upper_extents,
                                                 false);

      // merge the lower and upper mesh into one total_coarse_grid.
      // now we have at least two cells
      GridGenerator::merge_triangulations(lower_coarse_grid,
                                          upper_coarse_grid,
                                          total_coarse_grid);

      // set the boundary indicators
      set_boundary_indicators(total_coarse_grid);

      // tell p4est about the periodicity of the mesh.
      std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> >
      periodicity_vector;
      for (int i=0; i<dim+dim-1; ++i)
        {
          if (periodic[i])
            GridTools::collect_periodic_faces
            (total_coarse_grid, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
             /*direction*/ i%dim, periodicity_vector);
        }
      if (periodicity_vector.size() > 0)
        total_coarse_grid.add_periodicity (periodicity_vector);

      // make sure the right boundary indicators are set after refinement
      // through the function set_boundary_indicators above
      total_coarse_grid.signals.post_refinement.connect
      (std_cxx1x::bind (&TwoMergedBoxes<dim>::set_boundary_indicators,
                        std_cxx1x::cref(*this),
                        std_cxx1x::ref(total_coarse_grid)));
    }


    template <int dim>
    std::set<types::boundary_id>
    TwoMergedBoxes<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim+2*(dim-1)-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim+2*(dim-1); ++i)
        s.insert (i);
      return s;
    }



    template <int dim>
    std::map<std::string,types::boundary_id>
    TwoMergedBoxes<dim>::
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
                  std::pair<std::string,types::boundary_id>("top",    3),
                  std::pair<std::string,types::boundary_id>("left lithosphere", 4),
                  std::pair<std::string,types::boundary_id>("right lithosphere",5)
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
                  std::pair<std::string,types::boundary_id>("top",    5),
                  std::pair<std::string,types::boundary_id>("left lithosphere",  6),
                  std::pair<std::string,types::boundary_id>("right lithosphere", 7),
                  std::pair<std::string,types::boundary_id>("front lithosphere", 8),
                  std::pair<std::string,types::boundary_id>("back lithosphere",  9)
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
    TwoMergedBoxes<dim>::
    get_periodic_boundary_pairs () const
    {
      std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundaries;
      for ( unsigned int i=0; i<dim+dim-1; ++i)
        if (periodic[i])
          periodic_boundaries.insert( std::make_pair( std::pair<types::boundary_id, types::boundary_id>(2*i, 2*i+1), i) );
      return periodic_boundaries;
    }

    template <int dim>
    Point<dim>
    TwoMergedBoxes<dim>::get_extents () const
    {
      return extents;
    }

    template <int dim>
    Point<dim>
    TwoMergedBoxes<dim>::get_origin () const
    {
      return lower_box_origin;
    }

    template <int dim>
    double
    TwoMergedBoxes<dim>::
    length_scale () const
    {
      return 0.01*extents[0];
    }


    template <int dim>
    double
    TwoMergedBoxes<dim>::depth(const Point<dim> &position) const
    {
      const double d = maximal_depth()-(position(dim-1)-lower_box_origin[dim-1]);
      return std::min (std::max (d, 0.), maximal_depth());
    }


    template <int dim>
    Point<dim>
    TwoMergedBoxes<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // choose a point on the center axis of the domain
      Point<dim> p = extents/2+lower_box_origin;
      p[dim-1] = extents[dim-1]+lower_box_origin[dim-1]-depth;

      return p;
    }


    template <int dim>
    double
    TwoMergedBoxes<dim>::maximal_depth() const
    {
      return extents[dim-1];
    }

    template <int dim>
    bool
    TwoMergedBoxes<dim>::has_curved_elements() const
    {
      return false;
    }

    template <int dim>
    void
    TwoMergedBoxes<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          prm.declare_entry ("Lithospheric thickness", "0.2",
                             Patterns::Double (0),
                             "The thickness of the lithosphere used to create "
                             "additional boundary indicators to set specific "
                             "boundary conditions for the lithosphere. ");

          // Total box extents
          prm.declare_entry ("X extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in x-direction. Units: m.");
          prm.declare_entry ("Y extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in y-direction. Units: m.");
          prm.declare_entry ("Z extent", "1",
                             Patterns::Double (0),
                             "Extent of the box in z-direction. This value is ignored "
                             "if the simulation is in 2d. Units: m.");

          // Total box origin
          prm.declare_entry ("Box origin X coordinate", "0",
                             Patterns::Double (),
                             "X coordinate of box origin. Units: m.");
          prm.declare_entry ("Box origin Y coordinate", "0",
                             Patterns::Double (),
                             "Y coordinate of box origin. Units: m.");
          prm.declare_entry ("Box origin Z coordinate", "0",
                             Patterns::Double (),
                             "Z coordinate of box origin. This value is ignored "
                             "if the simulation is in 2d. Units: m.");

          // Lower box repetitions
          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction of the lower box. "
                             "The same number of repetitions will be used in the upper box.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction of the lower box. If the simulation "
                             "is in 3d, the same number of repetitions will be used in the upper box.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction of the lower box. "
                             "This value is ignored if the simulation is in 2d.");

          // Upper box repetitions
          prm.declare_entry ("Y repetitions lithosphere", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction in the lithosphere. "
                             "This value is ignored if the simulation is in 3d.");
          prm.declare_entry ("Z repetitions lithosphere", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction in the lithosphere. "
                             "This value is ignored if the simulation is in 2d.");

          // Whole box periodicity
          prm.declare_entry ("X periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in X direction.");
          prm.declare_entry ("Y periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Y direction.");
          prm.declare_entry ("Z periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Z direction. "
                             "This value is ignored if the simulation is in 2d.");
          prm.declare_entry ("X periodic lithosphere", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in X direction in the lithosphere.");
          prm.declare_entry ("Y periodic lithosphere", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Y direction in the lithosphere. "
                             "This value is ignored if the simulation is in 2d. ");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    TwoMergedBoxes<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box with lithosphere boundary indicators");
        {
          const double thickness_lith = prm.get_double("Lithospheric thickness");

          extents[0]           = prm.get_double ("X extent");
          lower_extents[0]     = extents[0];
          upper_extents[0]     = extents[0];
          lower_box_origin[0]  = prm.get_double ("Box origin X coordinate");
          upper_box_origin[0]  = lower_box_origin[0];
          periodic[0]          = prm.get_bool ("X periodic");
          periodic[dim]        = prm.get_bool ("X periodic lithosphere");
          // to match the two triangulations, it is required that the
          // number of horizontal repetitions are the same
          lower_repetitions[0] = prm.get_integer ("X repetitions");
          upper_repetitions[0] = lower_repetitions[0];

          extents[1]           = prm.get_double ("Y extent");
          lower_box_origin[1]  = prm.get_double ("Box origin Y coordinate");
          periodic[1]          = prm.get_bool ("Y periodic");;
          lower_repetitions[1] = prm.get_integer ("Y repetitions");

          if (dim == 2)
            {
              lower_extents[1]     = extents[1] - thickness_lith;
              upper_extents[1]     = thickness_lith;
              upper_box_origin[1]  = lower_extents[1];
              upper_repetitions[1] = prm.get_integer ("Y repetitions lithosphere");
            }

          if (dim == 3)
            {
              lower_extents[1]     = extents[1];
              upper_extents[1]     = extents[1];
              upper_box_origin[1]  = lower_box_origin[1];
              periodic[dim+1]      = prm.get_bool ("Y periodic lithosphere");
              // to match the two triangulations, it is required that the
              // number of horizontal repetitions are the same
              upper_repetitions[1] = lower_repetitions[1];
              extents[2]           = prm.get_double ("Z extent");
              lower_extents[2]     = extents[2] - thickness_lith;
              upper_extents[2]     = thickness_lith;
              lower_box_origin[2]  = prm.get_double ("Box origin Z coordinate");
              upper_box_origin[2]  = lower_extents[2];
              periodic[2]          = prm.get_bool ("Z periodic");
              lower_repetitions[2] = prm.get_integer ("Z repetitions");
              upper_repetitions[2] = prm.get_integer ("Z repetitions lithosphere");
            }

          height_lith = extents[dim-1] - thickness_lith;

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
    ASPECT_REGISTER_GEOMETRY_MODEL(TwoMergedBoxes,
                                   "box with lithosphere boundary indicators",
                                   "A box geometry parallel to the coordinate directions. "
                                   "The extent of the box in each coordinate direction "
                                   "is set in the parameter file. This geometry model labels its "
                                   "sides with 2*dim+2*(dim-1) boundary indicators: in 2d, boundary indicators 0 through 3 "
                                   "denote the left, right, bottom and top boundaries, while indicators"
                                   "4 and 5 denote the upper part of the left and right vertical boundary, "
                                   "respectively. In 3d, boundary "
                                   "indicators 0 through 5 indicate left, right, front, back, bottom "
                                   "and top boundaries (see also the documentation of the deal.II class "
                                   "``GeometryInfo''), while indicators 6, 7, 8 and 9 denote the left, "
                                   "rigth, front and back upper parts of the vertical boundaries, respectively. "
                                   "You can also use symbolic names ``left'', ``right'', "
                                   "``left lithosphere'', etc., to refer to these boundaries in input files."
                                   "\n\n"
                                   "Note that for a given ``Global refinement level'' and no user-specified "
                                   "``Repetitions'', the lithosphere part of the mesh will be more refined. "
                                   "\n\n"
                                   "The additional boundary indicators for the lithosphere allow for "
                                   "selecting boundary conditions for the "
                                   "lithosphere different from those for the underlying mantle. "
                                   "An example application of this geometry is to prescribe a velocity on "
                                   "the lithospheric plates, but use open boundary conditions underneath. ")
  }
}
