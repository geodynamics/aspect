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


#ifndef __aspect__geometry_model_ellipsoidal_chunk_h
#define __aspect__geometry_model_ellipsoidal_chunk_h

#include <aspect/geometry_model/interface.h>
#include <deal.II/grid/manifold.h>

/**
 * This geometry model implements an (3d) ellipsoidal chunk geometry where two of the axis have
 * the same length. The ellipsoidal chunk can be non-coordinate parallel part of the ellipsoid.
 * @author This plugin is a joined effort of Menno Fraters, D Sarah Stamps and Wolfgang Bangerth
 */

namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes a geometry for an ellipsoid such as the WGS84 model of the earth.
     */
    template <int dim>
    class EllipsoidalChunk : public Interface<dim>
    {
      public:
        /**
         * A class which describes the manifold.
         */
        class EllipsoidalChunkGeometry : public ChartManifold<dim,3,3>
        {
          public:

            void
            set_initial_values(double para_semi_major_axis_a,
                               double para_eccentricity,
                               double para_semi_minor_axis_b,
                               double para_rot_para_to_para_angle,
                               double para_para_to_rect_angle,
                               double para_rotation_longitude,
                               double para_rotation_latitude,
                               double para_bottom_depth);

            virtual
            Point<3>
            pull_back(const Point<3> &space_point) const;

            virtual
            Point<2>
            pull_back(const Point<2> &space_point) const;

            virtual
            Point<3>
            push_forward(const Point<3> &chart_point) const;

            Point<dim>
            rotate_rectangle_to_parallelogram(const Point<dim> &x) const;

            Point<dim>
            rotate_parallelogram(const Point<dim> &x) const;



          private:
            double semi_major_axis_a = -1;
            double eccentricity = -1;
            double semi_minor_axis_b = -1;
            double rot_para_to_para_angle = 0;
            double para_to_rect_angle = 0;
            double rotation_longitude = -1;
            double rotation_latitude = -1;
            double bottom_depth = -1;

            Point<3> push_forward_ellipsoid (const Point<3> &phi_theta_d, const double semi_major_axis_a, const double eccentricity) const;
            Point<3> pull_back_ellipsoid (const Point<3> &x, const double semi_major_axis_a, const double eccentricity) const;
        };


        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void
        create_coarse_mesh(parallel::distributed::Triangulation<dim> &coarse_grid) const;

        /**
         * Return the typical length scale one would expect of features in this geometry,
         * assuming realistic parameters.
         *
         * We return 1/20th of the distance from the two endpoints (vertex_0 and vertex_7) of the cell.
         */
        virtual
        double
        length_scale() const;

        virtual
        double
        depth(const Point<dim> &position) const;

        virtual Point<dim>
        representative_point(const double depth) const;

        virtual
        double
        maximal_depth() const;

        /**
         * Return the set of boundary indicators that are used by this model. This
         * information is used to determine what boundary indicators can be used in
         * the input file.
         *
         * The box model uses boundary indicators zero through 2*dim-1, with the first
         * two being the faces perpendicular to the x-axis, the next two perpendicular
         * to the y-axis, etc.
         */
        virtual std::set<types::boundary_id>
        get_used_boundary_indicators() const;

        /*
        *Set symbolic names for boudaries (mrtf)
        */
        virtual std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map() const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters(ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters(ParameterHandler &prm);

        /**
         * Calculate radius at current position.
         */
        double
        get_radius(const Point<dim> &point) const;

        /**
         * Retrieve the semi minor axis b value.
         */
        double
        get_semi_minor_axis_b() const;

        /**
         * Retrieve the semi major axis a value.
         */
        double
        get_semi_major_axis_a() const;


        /**
         * Retrieve the value of the eccentricity.
         */
        double
        get_eccentricity() const;


        /**
         * Retrieve the rot_para_to_para_angle variable.
         */
        double
        get_rot_para_to_para_angle() const;


        /**
         * Retrieve the para_to_rect_angle variable.
         */
        double
        get_para_to_rect_angle() const;


        /**
         * Retrieve the rotation_longitude variable.
         */
        double
        get_rotation_longitude() const;


        /**
         * Retrieve the rotation_latitude variable.
         */
        double
        get_rotation_latitude() const;


        /**
         * Retrieve the manifold object.
         */
        EllipsoidalChunkGeometry
        get_manifold() const;

      private:
        // Declare variables for reading in coordinates of the region of interest.
        std::vector<std::vector<double> > corners;
        std::vector<std::vector<double> > corners_parallelogram;
        std::vector<std::vector<double> > corners_rectangle;
        double semi_major_axis_a;
        double eccentricity;
        double semi_minor_axis_b;
        double rot_para_to_para_angle;
        double para_to_rect_angle;
        double rotation_longitude;
        double rotation_latitude;
        double bottom_depth;
        double westLongitude;
        double eastLongitude;
        double northLatitude;
        double southLatitude;
        // Declare variables for subdividing
        int EW_subdiv;
        int NS_subdiv;
        int depth_subdiv;



        /**
         * Construct manifold object Pointer to an object that describes the geometry.
         */
        EllipsoidalChunkGeometry   manifold;

        static void set_manifold_ids (Triangulation<dim> &triangulation)
        {
          for (typename Triangulation<dim>::active_cell_iterator cell =
                 triangulation.begin_active(); cell != triangulation.end(); ++cell)
            cell->set_all_manifold_ids (15);
        }


        static void clear_manifold_ids (Triangulation<dim> &triangulation)
        {
          for (typename Triangulation<dim>::active_cell_iterator cell =
                 triangulation.begin_active(); cell != triangulation.end(); ++cell)
            cell->set_all_manifold_ids (numbers::invalid_manifold_id);
        }
    };
  }
}


#endif
