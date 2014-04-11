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
/*  $Id: box.cc 1538 2013-01-06 03:12:23Z bangerth $  */


#include <aspect/geometry_model/rebound.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/polynomial.h>
#include <cmath>
#include <vector>

namespace aspect
{
  template <int dim, int spacedim> class PerturbedBoundary;
}

namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    ReboundBox<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
     std::vector<unsigned int> repetitions(dim);
      //find minimum extent
      double min =extents[0];
      for(int i = 1; i< dim; ++i) if (extents[i] < min) min = extents[i];
      //get repetitions
      for(int i=0; i<dim; ++i) repetitions[i] = (unsigned int)(multiplicity*extents[i]/min);
      

      //make the regular grid
      GridGenerator::subdivided_hyper_rectangle (coarse_grid, repetitions,
                                      Point<dim>(),
                                      extents, true);
  

      //move the vertices
      std::vector<bool> vertex_touched (coarse_grid.n_vertices(), false);

      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
      for (cell = coarse_grid.begin_active();  cell != coarse_grid.end();  ++cell)
        for(unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;  ++v)
          if(vertex_touched[cell->vertex_index(v)] == false)
          {
            Point<dim> &vertex = cell->vertex(v);
            vertex[dim-1] = vertex[dim-1] + std::cos(order*2.0*M_PI*vertex[0]/extents[0])*
                                            amplitude*vertex[dim-1]/extents[dim-1];
            vertex_touched[cell->vertex_index(v)] = true;
          }

/*      GridGenerator::hyper_rectangle (coarse_grid,
                                      Point<dim>(),
                                      extents);
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        coarse_grid.begin_active()->face(f)->set_boundary_indicator(f);

      //Tell p4est about the periodicity of the mesh.  We make all directions but
      //the up down direction periodic in space
      std::vector<std_cxx1x::tuple< typename parallel::distributed::Triangulation<dim>::cell_iterator, unsigned int,
                                    typename parallel::distributed::Triangulation<dim>::cell_iterator, unsigned int> >
                                   periodicity_vector;
      for( unsigned int i=0; i<dim-1; ++i)
        GridTools::identify_periodic_face_pairs(coarse_grid, 2*i, 2*i+1, i, periodicity_vector);

      coarse_grid.add_periodicity(periodicity_vector);

      static const PerturbedBoundary<dim> bound(extents, order, amplitude);
      coarse_grid.set_boundary((dim == 2 ? 3 : 5) , bound);*/
    }


    template <int dim>
    std::set<types::boundary_id>
    ReboundBox<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }
/*
    template <int dim>
    std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
    ReboundBox<dim>::
    get_periodic_boundary_pairs () const
    {
      std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> > periodic_boundaries;
      for( unsigned int i=0; i<dim-1; ++i)
        periodic_boundaries.insert( std::make_pair( std::pair<types::boundary_id, types::boundary_id>(2*i, 2*i+1), i) );
      return periodic_boundaries;
    }
*/

    template <int dim>
    Point<dim>
    ReboundBox<dim>::get_extents () const
    {
      return extents;
    }


    template <int dim>
    double
    ReboundBox<dim>::
    length_scale () const
    {
      return 0.01*extents[0];
    }


    template <int dim>
    double
    ReboundBox<dim>::depth(const Point<dim> &position) const
    {
      const double d = maximal_depth()-position(dim-1);

      return d;
    }


    template <int dim>
    Point<dim>
    ReboundBox<dim>::representative_point(const double depth) const
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
    ReboundBox<dim>::maximal_depth() const
    {
      return extents[dim-1]*(2.0);
    }

    template <int dim>
    double
    ReboundBox<dim>::get_wavelength () const
    {
      return extents[0]/order;
    }

    template <int dim>
    void
    ReboundBox<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Rebound Box");
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
          prm.declare_entry ("Order", "1",
                             Patterns::Double (0),
                             "Order of the perturbation");
          prm.declare_entry ("Number of cells", "20",
                             Patterns::Integer (0),
                             "Number of cells in the smallest direction");
          prm.declare_entry ("Amplitude", "0.1",
                             Patterns::Double (0),
                             "Scaled amplitude of the perturbation");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ReboundBox<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Rebound Box");
        {
          extents[0] = prm.get_double ("X extent");

          if (dim >= 2)
            extents[1] = prm.get_double ("Y extent");

          if (dim >= 3)
            extents[2] = prm.get_double ("Z extent");
          order = prm.get_double ("Order");
          multiplicity = prm.get_integer ("Number of cells");
          amplitude = prm.get_double ("Amplitude");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ReboundSphere<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      if(R0 == 0.0)
      {
        GridGenerator::hyper_ball(coarse_grid, Point<dim>(), R1);
        static const HyperShellBoundary<dim> boundary_shell;
        coarse_grid.set_boundary (0, boundary_shell);
        coarse_grid.refine_global(refines);
        coarse_grid.set_boundary(0);
      }
      else
      {
        GridGenerator::hyper_shell (coarse_grid,
                                   Point<dim>(),
                                   R0,
                                   R1,
                                   (dim==3) ? 96 : 12,
                                   true);

        static const HyperShellBoundary<dim> boundary_shell;
        coarse_grid.set_boundary (0, boundary_shell);
        coarse_grid.set_boundary (1, boundary_shell);
        coarse_grid.refine_global(refines);
        coarse_grid.set_boundary(0);
        coarse_grid.set_boundary(1);
      }
 
      typename Triangulation<dim>::active_cell_iterator cell;
      std::vector< CellData<dim> > cells(0);
      std::vector<bool> vertex_touched (coarse_grid.n_vertices(), false);
        
      Polynomials::Legendre legendre(degree);

      for (cell = coarse_grid.begin_active(); cell != coarse_grid.end();  ++cell)
      {
        CellData<dim> c;
        for (unsigned int v = 0; v != GeometryInfo<dim>::vertices_per_cell;  ++v)
          if(vertex_touched[ cell->vertex_index(v)] == false)
          {
            c.vertices[v] = cell->vertex_index(v);
            Point<dim> &vertex = cell->vertex(v);
            double r = vertex.norm();
            double phi = std::atan2(vertex[1],vertex[0]);
            if(dim == 3) 
            {
              double theta = std::acos(vertex[2]/r); 
              double lval = legendre.value((std::cos(theta)+1.0)/2.0)/legendre.value(0.0);
              vertex = vertex + amplitude * ( (r-R0)/(R1-R0) ) * lval * vertex/vertex.norm();
            }
            else vertex = vertex + amplitude * ( (r-R0)/(R1-R0) ) * std::cos(degree*phi) * vertex/vertex.norm();
            vertex_touched[ cell->vertex_index(v)] = true;
          }
        cells.push_back(c);
      }
      std::vector<Point<dim> > vertices = coarse_grid.get_vertices();
    }


    template <int dim>
    std::set<types::boundary_id>
    ReboundSphere<dim>::
    get_used_boundary_indicators () const
    {
      const types::boundary_id s[] = { 0, 1 };
      return std::set<types::boundary_id>(&s[0],
                                            &s[sizeof(s)/sizeof(s[0])]);
    }


    template <int dim>
    double
    ReboundSphere<dim>::
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
    ReboundSphere<dim>::depth(const Point<dim> &position) const
    {
      return std::min (std::max (R1-position.norm(), 0.), maximal_depth());
    }

    template <int dim>
    unsigned int
    ReboundSphere<dim>::get_degree( ) const
    {
      return degree;
    }


    template <int dim>
    Point<dim>
    ReboundSphere<dim>::representative_point(const double depth) const
    {
      Point<dim> p;
      p(dim-1) = std::min (std::max(R1 - depth, R0), R1);
      return p;
    }



    template <int dim>
    double
    ReboundSphere<dim>::maximal_depth() const
    {
      return (R1-R0)*1.5;
    }

    template <int dim>
    double ReboundSphere<dim>::inner_radius () const
    {
      return R0;
    }



    template <int dim>
    double ReboundSphere<dim>::outer_radius () const
    {
      return R1;
    }

    template <int dim>
    void
    ReboundSphere<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Rebound sphere");
        {
          prm.declare_entry ("Inner radius", "3481000",  // 6371-2890 in km
                             Patterns::Double (0),
                             "Inner radius of the spherical shell. Units: m.");
          prm.declare_entry ("Outer radius", "6336000",  // 6371-35 in km
                             Patterns::Double (0),
                             "Outer radius of the spherical shell. Units: m.");
          prm.declare_entry ("Degree", "2",  // 6371-35 in km
                             Patterns::Integer (0),
                             "Order of perterbation.");
          prm.declare_entry ("Refines", "4",  // 6371-35 in km
                             Patterns::Integer (0),
                             "Order of perterbation.");
          prm.declare_entry ("Amplitude", "100000",  // 6371-35 in km
                             Patterns::Double (0),
                             "Amplitude of perturbation. Units: m.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ReboundSphere<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Rebound sphere");
        {
          R0            = prm.get_double ("Inner radius");
          R1            = prm.get_double ("Outer radius");
          degree = prm.get_integer("Degree");
          refines = prm.get_integer("Refines");
          amplitude     = prm.get_double ("Amplitude");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

namespace aspect
{
   using namespace dealii; 
 
  template <int dim, int spacedim=dim>
  class PerturbedBoundary : public StraightBoundary<dim, spacedim>
  {
    public:
      PerturbedBoundary (const Point<spacedim> grid_extents,
                         const unsigned int perturbation_order,
                         const double perturbation_amplitude):
                         extents(grid_extents),
                         order(perturbation_order),
                         amplitude(perturbation_amplitude) {}
                         

      virtual Point<spacedim>
      get_new_point_on_line (const typename Triangulation<dim,spacedim>::line_iterator &line) const
      {
        Point<spacedim> middle = StraightBoundary<dim,spacedim>::get_new_point_on_line (line);
 
        if(dim == 2)
          middle[1] += amplitude*std::sin(2.0*M_PI*order*middle[0]/extents[0]);
        else if (dim == 3)
          middle[2] += amplitude*std::sin(2.0*M_PI*order*middle[0]/extents[0])*
                                std::sin(2.0*M_PI*order*middle[1]/extents[1]);
        return middle;
      }

      virtual Point<spacedim>
      get_new_point_on_quad (const typename Triangulation<dim,spacedim>::quad_iterator &quad) const
      {
        Point<dim> middle = StraightBoundary<dim,spacedim>::get_new_point_on_quad (quad);
        if(dim == 2)
          middle[1] += amplitude*std::sin(2.0*M_PI*order*middle[0]/extents[0]);
        else if (dim == 3)
          middle[2] += amplitude*std::sin(2.0*M_PI*order*middle[0]/extents[0])*
                                std::sin(2.0*M_PI*order*middle[1]/extents[1]);
        return middle;
      }

/*      virtual void
      get_intermediate_points_on_line (const typename Triangulation<dim>::line_iterator &line,
                                       std::vector<Point<dim> > &points) const;

      virtual void
      get_intermediate_points_on_quad (const typename Triangulation<dim>::quad_iterator &quad,
                                       std::vector<Point<dim> > &points) const;
*/
      virtual void
      get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &face,
                               typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const
      {
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
        {
          const Point<spacedim> vertex = face->vertex(v);
          Point<spacedim> normal;
          if(dim == 2)
          {
            normal = Point<dim>::unit_vector(1);
            normal[0] -= amplitude*amplitude*order*2.0*M_PI/extents[0]
                         *std::sin(2.0*M_PI*order*vertex[0]/extents[0]);
          }
          else if(dim == 3)
          {
            normal = Point<dim>::unit_vector(2);
            normal[0] -= amplitude*amplitude*order*2.0*M_PI/extents[0]
                         *std::sin(2.0*M_PI*order*vertex[0]/extents[0])
                         *std::sin(2.0*M_PI*order*vertex[1]/extents[1]);
            normal[1] -= amplitude*amplitude*order*2.0*M_PI/extents[1]
                         *std::sin(2.0*M_PI*order*vertex[0]/extents[0])
                         *std::sin(2.0*M_PI*order*vertex[1]/extents[1]);
          }
          face_vertex_normals[v] = (normal/normal.norm());
        }
      } 
    private:
      const Point<spacedim> extents;
      const double amplitude;
      const unsigned int order;
  };
}

// explicit instantiations
namespace aspect
{
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(ReboundBox,
                                   "rebound box",
                                   "Geometry model for benchmarking isostatic rebound.")
    ASPECT_REGISTER_GEOMETRY_MODEL(ReboundSphere,
                                   "rebound sphere",
                                   "")
  }
}
