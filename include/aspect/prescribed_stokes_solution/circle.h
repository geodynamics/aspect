// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
