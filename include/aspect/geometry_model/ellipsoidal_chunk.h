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


#ifndef _aspect_geometry_model_ellipsoidal_chunk_h
#define _aspect_geometry_model_ellipsoidal_chunk_h

#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/grid/manifold.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes a geometry for an ellipsoid such as the WGS84 model of the earth.
     *
     * This geometry model implements a (3d) ellipsoidal chunk geometry where two of the axis have
     * the same length. The ellipsoidal chunk can be a non-coordinate parallel part of the ellipsoid.
     *
     * @author This plugin is a joined effort of Menno Fraters, D. Sarah Stamps and Wolfgang Bangerth
     */
    template <int dim>
    class EllipsoidalChunk : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * A class which describes the manifold.
         */
        class EllipsoidalChunkGeometry : public ChartManifold<dim,3,3>
        {
          public:
            /**
             * Constructor
             */
            EllipsoidalChunkGeometry();

            /**
             * Copy constructor
             */
            EllipsoidalChunkGeometry(const EllipsoidalChunkGeometry &other);

            /**
             * An initialization function necessary to make sure that the
             * manifold has access to the topography plugins.
             */
            void
            initialize(const InitialTopographyModel::Interface<dim> *topography);

            /**
             * Sets several parameters for the ellipsoidal manifold object.
             */
            void
            set_manifold_parameters(const double para_semi_major_axis_a,
                                    const double para_eccentricity,
                                    const double para_semi_minor_axis_b,
                                    const double para_bottom_depth,
                                    const std::vector<Point<2> > &para_corners);

            /**
             * The deal.ii pull back function in 3d. This function receives
             * cartesian points x,y and z and returns spherical/ellipsoidal
             * coordinates phi, theta and depth, also accounting for the
             * topography.
             */
            Point<3>
            pull_back(const Point<3> &space_point) const override;

            /**
             * The deal.ii pull back function in 2d. This function should
             * not be used, until the TODO in the cc file has been fixed.
             */
            virtual
            Point<2>
            pull_back(const Point<2> &space_point) const;

            /**
             * The deal.ii push forward function in 3d. This function receives
             * spherical/ellipsoidal coordinates phi, theta and depth and
             * returns cartesian points x,y and z, also accounting for the
             * topography.
             */
            Point<3>
            push_forward(const Point<3> &chart_point) const override;

            /**
             * This function does the actual pull back from the ellipsoid.
             * For the equation details, please see deal.ii step 53.
             */
            Point<3> pull_back_ellipsoid (const Point<3> &x, const double semi_major_axis_a, const double eccentricity) const;

            /**
             * This function does the actual push forward to the ellipsoid.
             * For the equation details, please see deal.ii step 53.
             */
            Point<3> push_forward_ellipsoid (const Point<3> &phi_theta_d, const double semi_major_axis_a, const double eccentricity) const;

            /**
             * Return a copy of this manifold.
             */
            std::unique_ptr<Manifold<dim,3> >
            clone() const override;

          private:
            /**
             * This function adds topography to the cartesian coordinates.
             * For the equation details, please see deal.ii step 53.
             */
            Point<3> push_forward_topography (const Point<3> &phi_theta_d_hat) const;

            /**
             * This function removes topography from the cartesian coordinates.
             * For the equation details, please see deal.ii step 53.
             */
            Point<3> pull_back_topography (const Point<3> &phi_theta_d) const;


            double semi_major_axis_a;
            double eccentricity;
            double semi_minor_axis_b;
            double bottom_depth;
            std::vector<Point<2> > corners;
            const InitialTopographyModel::Interface<dim> *topography;
        };

        /**
         * Initialize function
         */
        void
        initialize () override;


        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        void
        create_coarse_mesh(parallel::distributed::Triangulation<dim> &coarse_grid) const override;

        /**
         * Return the typical length scale one would expect of features in this geometry,
         * assuming realistic parameters.
         */
        double
        length_scale() const override;

        /**
         * Return the depth that corresponds to the given
         * position. The documentation of the base class (see
         * GeometryModel::Interface::depth()) describes in detail how
         * "depth" is interpreted in general.
         *
         * Computing a depth requires a geometry model to define a
         * "vertical" direction. The current class considers the
         * radial vector away from the origin as vertical and
         * considers the "outer" boundary as the "surface". In almost
         * all cases one will use a gravity model that also matches
         * these definitions.
         */
        double
        depth(const Point<dim> &position) const override;

        /**
         * Placeholder for a function returning the height of the given
         * position relative to the reference model surface.
         */
        double
        height_above_reference_surface(const Point<dim> &position) const override;

        /**
         * Returns a point in the center of the domain.
         */
        Point<dim>
        representative_point(const double depth) const override;

        /**
         * Return whether the given point lies within the domain specified
         * by the geometry. This function does not take into account
         * initial or dynamic surface topography.
         */
        bool
        point_is_in_domain(const Point<dim> &point) const override;

        /**
         * Returns the bottom depth which was used to create the geometry and
         * which is defined by the depth parameter.
         */
        double
        maximal_depth() const override;

        /**
         * Return the set of boundary indicators that are used by this model. This
         * information is used to determine what boundary indicators can be used in
         * the input file.
         *
         * The box model uses boundary indicators zero through 2*dim-1, with the first
         * two being the faces perpendicular to the x-axis, the next two perpendicular
         * to the y-axis, etc.
         */
        std::set<types::boundary_id>
        get_used_boundary_indicators() const override;

        /**
         * Set symbolic names for boundaries (mrtf)
         */
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map() const override;

        /**
         * Returns what the natural coordinate system for this geometry model is,
         * which for a Ellipsoidal chunk is Ellipsoidal.
         */
        aspect::Utilities::Coordinates::CoordinateSystem natural_coordinate_system() const override;

        /**
         * Takes the Cartesian points (x,z or x,y,z) and returns standardized
         * coordinates which are most 'natural' to the geometry model. For a
         * ellipsoidal chunk this is (radius, longitude) in 2d and (radius,
         * longitude, latitude) in 3d. Note that internally the coordinates are
         * stored in longitude, latitude, depth.
         */
        std::array<double,dim> cartesian_to_natural_coordinates(const Point<dim> &position) const override;

        /**
         * Undoes the action of cartesian_to_natural_coordinates, and turns the
         * coordinate system which is most 'natural' to the geometry model into
         * Cartesian coordinates.
         */
        Point<dim> natural_to_cartesian_coordinates(const std::array<double,dim> &position) const override;

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
        void
        parse_parameters(ParameterHandler &prm) override;

        /**
         * Calculate radius at current position.
         */
        double
        get_radius(const Point<dim> &position) const;

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
         * Retrieve the corners used to create the ellipsoid. This variable
         * contains four vectors with each two elements. Each set of two
         * elements represents a longitude and latitude value. The four
         * vectors represent respectively the point in the North-East,
         * North-West, South-West and South-East.
         */
        const std::vector<Point<2> > &
        get_corners() const;


        /**
         * Retrieve the manifold object.
         */
        EllipsoidalChunkGeometry
        get_manifold() const;

      private:
        // Declare variables for reading in coordinates of the region of interest.
        std::vector<Point<2> > corners;
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
        unsigned int EW_subdiv;
        unsigned int NS_subdiv;
        unsigned int depth_subdiv;



        /**
         * Construct manifold object Pointer to an object that describes the geometry.
         */
        EllipsoidalChunkGeometry   manifold;

        void
        set_boundary_ids(parallel::distributed::Triangulation<dim> &coarse_grid) const;
    };
  }
}


#endif
