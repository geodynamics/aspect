/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_lithosphere_rift_h
#define _aspect_initial_composition_lithosphere_rift_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the composition
     * considering a lithosphere consisting of an upper crust,
     * lower crust and mantle part.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class LithosphereRift : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const;

        /**
         * Return the overall shortest distance to the rift segments
         */
        double distance_to_rift (const Point<dim-1> &position) const;

        /*
         * Return the overall shortest distance to the polygon segments
         */
        std::pair<double,unsigned int> distance_to_polygon (const Point<dim-1> &position) const;

        /*
         * Return the position of the point in surface coordinates
         * i.e. x(,y) in meters or lon(,lat) in degrees.
         */
        Point<dim-1> surface_position (const Point<dim> &position,
                                       const bool cartesian_geometry) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * The standard deviation of the Gaussian amplitude of the lithospheric thicknesses
         * with distance from the rift axis.
         */
        double sigma_rift;

        /**
         * The half width of the hyperbolic tangent used to smooth the transitions
         * between reference and polygon lithospheric thicknesses.
         */
        double sigma_polygon;

        /**
         * Whether or not to take the polygon thicknesses as dominant, or to smooth them
         * gradually into rift areas.
         */
        bool blend_rift_and_polygon;

        /**
         * The maximum amplitude of the Gaussian distribution of the lithospheric thicknesses
         * with distance from the rift axis. It should have values between -1 and 1, where positive
         * numbers represent a reduction in thickness and negative numbers an increase.
         */
        std::vector<double> A;

        /**
         * The list of line segments consisting of two 2d coordinates per segment.
         * The segments represent the rift axis.
         */
        std::vector<std::array<Point<2>,2 > > point_list;

        /**
         * The list of lists of polygon points.
         * The polygon represents an area of different lithospheric thicknesses.
         */
        std::vector<std::vector<Point<2> > > polygon_point_list;

        /**
         * Vector for the reference field thicknesses away from the rift.
         */
        std::vector<double> thicknesses;

        /**
         * Vector for the field thicknesses inside the polygons.
         */
        std::vector<std::vector<double> > polygon_thicknesses;

    };
  }
}


#endif
