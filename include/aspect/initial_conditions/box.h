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
/*  $Id$  */


#ifndef __aspect__initial_conditions_box_h
#define __aspect__initial_conditions_box_h

#include <aspect/initial_conditions/interface.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes a perturbed initial temperature field for
     * a box geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class PerturbedBox : public Interface<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;
    };
  }
}

#endif
