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
#include <aspect/geometry_model/box.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria.h>
#include <cmath>
#include <vector>

namespace aspect
{
  namespace GeometryModel
  {
    /**
     * A class deriving from Box<dim>, which changes the upper boundary 
       with a sinusoidal perturbation of a given order and amplitude
     */
    template <int dim>
    class ReboundBox : public Box<dim>
    {
      public:
        /**
         * Generate a coarse mesh for the geometry described by this class.
         * Makes perturbs the top boundary of the box with a function
         * of the form z' = amplitude * cos(order * x )
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Give the depth of a point.
         */
        virtual
        double depth( const Point<dim> &position) const;

        /**
         * Give the maximal depth of a point.
         */
        virtual
        double maximal_depth() const;

      private:

        unsigned int order;  //Order of the perturbation
        double amplitude;  //amplitude of the perturbation

    };

    template <int dim>
    void
    ReboundBox<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
  
      //Call the normal Box mesh generator
      Box<dim>::create_coarse_mesh( coarse_grid );

      //move the vertices
      std::vector<bool> vertex_touched (coarse_grid.n_vertices(), false);

      typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell;
      for (cell = coarse_grid.begin_active();  cell != coarse_grid.end();  ++cell)
        for(unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;  ++v)
          if(vertex_touched[cell->vertex_index(v)] == false)
          {
            Point<dim> &vertex = cell->vertex(v);
            vertex[dim-1] = vertex[dim-1] + std::cos(order*2.0*M_PI*vertex[0]/(this->get_extents()[0]))*
                                            amplitude*vertex[dim-1]/(this->get_extents()[dim-1]);
            vertex_touched[cell->vertex_index(v)] = true;
          }

    }

    template <int dim>
    double
    ReboundBox<dim>::maximal_depth() const
    {
      return Box<dim>::maximal_depth()+amplitude;
    }

    template <int dim>
    double
    ReboundBox<dim>::depth(const Point<dim> &position) const
    {
      double d = maximal_depth()-position(dim-1);
      return d;
    }

    template <int dim>
    void
    ReboundBox<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      Box<dim>::declare_parameters(prm);

      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Rebound Box");
        {
          prm.declare_entry ("Order", "1",
                             Patterns::Double (0),
                             "Order of the perturbation");
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
      Box<dim>::parse_parameters(prm);

      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Rebound Box");
        {
          order = prm.get_double ("Order");
          amplitude = prm.get_double ("Amplitude");
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
    ASPECT_REGISTER_GEOMETRY_MODEL(ReboundBox,
                                   "rebound box",
                                   "Geometry model for benchmarking relaxation of topography.")
  }
}
