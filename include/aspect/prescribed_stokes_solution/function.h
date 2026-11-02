// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_prescribed_stokes_solution_function_h
#define _aspect_prescribed_stokes_solution_function_h

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    /**
     * A class that implements velocity and pressure solutions based on a
     * functional description provided in the input file.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function ();

        /**
         * Return the velocity and pressure as a function of position.
         */
        void
        stokes_solution (const Point<dim> &position, Vector<double> &value) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the velocity and
         * pressure values will next be evaluated. For the current class,
         * the function passes to the parsed function what the current time is.
         */
        void
        update () override;

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
         * A function object representing the components of the velocity.
         */
        Functions::ParsedFunction<dim> prescribed_velocity_function;
        /**
         * A function object representing the pressure.
         */
        Functions::ParsedFunction<dim> prescribed_pressure_function;
        /**
         * A function object representing the fluid pressure (in models with melt transport).
         */
        Functions::ParsedFunction<dim> prescribed_fluid_pressure_function;
        /**
         * A function object representing the compaction pressure (in models with melt transport).
         */
        Functions::ParsedFunction<dim> prescribed_compaction_pressure_function;
        /**
         * A function object representing the components of the fluid velocity (in models with melt transport).
         */
        Functions::ParsedFunction<dim> prescribed_fluid_velocity_function;
        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
