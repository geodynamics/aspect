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


#ifndef __aspect__boundary_composition_two_merged_boxes_h
#define __aspect__boundary_composition_two_merged_boxes_h

#include <aspect/boundary_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a
     * box geometry with additional boundary indicators for the
     * lithospheric part of the left and right (and front and back in 3D) boundaries.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class TwoMergedBoxes : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns constant compositions at the boundaries.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        virtual
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &location,
                                     const unsigned int compositional_field) const;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary compositions.
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

        /**
         * This function performs some basic sanity checks on the parameter
         * values previously read from the input file.
         */
        virtual void initialize ();

      private:
        /**
         * The values of the various composition variables on each of the
         * 2*dim+2*(dim-1) boundaries of the box.
         */
        std::vector<double> composition_values[2*dim+2*(dim-1)];
    };
  }
}


#endif
