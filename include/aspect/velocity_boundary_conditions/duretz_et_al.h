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


#ifndef __aspect__velocity_boundary_conditions_duretz_et_al_h
#define __aspect__velocity_boundary_conditions_duretz_et_al_h

#include <aspect/velocity_boundary_conditions/interface.h>

namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    namespace DuretzEtAl
    {
      using namespace dealii;

      /**
       * A class that implements the velocity boundary conditions for the
       * ``inclusion'' benchmark. See the manual for more information on this
       * benchmark.
       *
       * @ingroup VelocityBoundaryConditionsModels
       */
      template <int dim>
      class Inclusion : public Interface<dim>
      {
        public:
          /**
           * Constructor.
           */
          Inclusion();

          /**
           * Return the boundary velocity as a function of position.
           */
          virtual
          Tensor<1,dim>
          boundary_velocity (const Point<dim> &position) const;

        private:
          double eta_B;
      };
    }
  }
}


#endif
