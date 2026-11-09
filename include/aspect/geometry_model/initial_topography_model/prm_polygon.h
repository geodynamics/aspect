// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_geometry_model_initial_topography_prm_polygon_h
#define _aspect_geometry_model_initial_topography_prm_polygon_h

#include <aspect/geometry_model/initial_topography_model/interface.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
     * A class that describes an initial topography for the geometry model,
     * by defining a set of polygons on the surface from the prm file. It
     * sets the elevation in each polygon to a constant value.
     */
    template <int dim>
    class PrmPolygon : public Interface<dim>
    {
      public:
        /**
         * Return the value of the topography for a point.
         */
        double value (const Point<dim-1> &p) const override;

        /**
         * Return the maximum value of the elevation.
         */
        double max_topography () const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The values of the topography are stored in a vector.
         */
        std::vector<double> topography_values;

        /**
         * The maximum topography in this model
         */
        double maximum_topography;

        /**
         * The polygons and their points are stored in this vector.
         */
        std::vector<std::vector<Point<2>>> point_lists;

    };
  }
}


#endif
