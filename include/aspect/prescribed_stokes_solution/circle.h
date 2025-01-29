/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_prescribed_stokes_solution_circle_h
#define _aspect_prescribed_stokes_solution_circle_h

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace PrescribedStokesSolution
  {
    /**
     * A class that implements a circular, divergence-free flow field
     * around the origin of the coordinate system.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    class Circle : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        void stokes_solution (const Point<dim> &p, Vector<double> &value) const override;
    };
  }
}


#endif
