/*
  Copyright (C) 2013 by the authors of the ASPECT code.

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
/*  $Id: box.h 1538 2013-01-06 03:12:23Z bangerth $  */


#ifndef __aspect__boundary_composition_box_h
#define __aspect__boundary_composition_box_h

#include <aspect/boundary_composition/interface.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a box
     * geometry.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class Box : public Interface<dim>
    {
      public:
        /**
         * Return the composition that is to hold at a particular location on the
         * boundary of the domain. This function returns constant compositions
         * at the left and right boundaries.
         *
         * @param geometry_model The geometry model that describes the domain. This may
         *   be used to determine whether the boundary composition model is implemented
         *   for this geometry.
         * @param boundary_indicator The boundary indicator of the part of the boundary
         *   of the domain on which the point is located at which we are requesting the
         *   composition.
         * @param location The location of the point at which we ask for the composition.
         **/
        virtual
        double composition (const GeometryModel::Interface<dim> &geometry_model,
                            const unsigned int                   boundary_indicator,
                            const Point<dim>                    &location,
                            const unsigned int                   compositional_field) const;

        /**
         * Declare the parameters this class takes through input files.
         * This class declares the inner and outer boundary compositions.
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

      private:
        double composition_[2*dim];
    };
  }
}


#endif
