/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model__initial_topography_prm_polygon_h
#define _aspect_geometry_model__initial_topography_prm_polygon_h

#include <aspect/geometry_model/initial_topography_model/interface.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    using namespace dealii;

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
        virtual
        double value (const Point<dim-1> &p) const;


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
         * The values of the topography are stored in a vector.
         */
        std::vector<double> topography_values;

        /**
         * The polygons and their points are stored in this vector.
         */
        std::vector<std::vector<Point<2> > > point_lists;

    };
  }
}


#endif
