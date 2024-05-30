/*
  Copyright (C) 2023 by the authors of the ASPECT code.

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
        initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const override;

        /**
         * Return the overall shortest distance to the rift center segments.
         */
        double distance_to_rift (const Point<dim-1> &position) const;

        /*
         * Return the overall shortest distance to the polygon segments. Polygons
         * can be used to specify lithosphere of different thicknesses.
         */
        std::pair<double,unsigned int> distance_to_polygon (const Point<dim-1> &position) const;

        /*
         * Return the position of the point in surface coordinates,
         * i.e. x(,y) in meters or lon(,lat) in degrees.
         */
        Point<dim-1> surface_position (const Point<dim> &position,
                                       const bool cartesian_geometry) const;

        /*
         * Compute the local thicknesses of the upper and lower crust
         * and lithospheric mantle based on the distance to the rift center
         * and the edge of polygons.
         */
        std::vector<double> compute_local_thicknesses(const Point<dim-1> &position) const;

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
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * The standard deviation of the rift-perpendicular Gaussian distribution
         * of the thinning/thickening of the lithospheric thicknesses with its
         * maximum along the rift axis.
         */
        double sigma_rift;

        /**
         * The maximum amplitude of the Gaussian distribution of the thinning/thickening
         * of the lithospheric thicknesses with distance from the rift axis.
         * The amplitude should have values between -1 and 1, where positive
         * numbers represent a reduction in thickness and negative numbers an increase.
         * For example, values of 0.25, 0, 0 reduce the reference thickness of the
         * upper crust by 25%, while the lower crust and mantle lithosphere are
         * untouched.
         */
        std::vector<double> A_rift;

        /**
         * The half width of the hyperbolic tangent used to smooth the transitions
         * between reference and polygon lithospheric thicknesses.
         */
        double sigma_polygon;

        /**
         * The list of line segments consisting of two 2d coordinates per segment.
         * The segments represent the rift axis.
         */
        std::vector<std::array<Point<2>,2 >> rift_point_list;

        /**
         * The list of lists of polygon points.
         * The polygons can represent areas of different lithospheric thicknesses.
         */
        std::vector<std::vector<Point<2>>> polygon_point_list;

        /**
         * Vector for the reference field thicknesses away from rift and polygons.
         */
        std::vector<double> reference_thicknesses;

        /**
         * Vector for the field thicknesses for each polygon.
         */
        std::vector<std::vector<double>> polygon_thicknesses;
    };
  }
}


#endif
