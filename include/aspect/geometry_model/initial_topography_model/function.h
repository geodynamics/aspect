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


#ifndef _aspect_geometry_model_initial_topography_model_function_h
#define _aspect_geometry_model_initial_topography_model_function_h

#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
     * A class that implements initial topography based
     * on a user-defined function..
     *
     * @ingroup InitialTopographyModels
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
         * Return the value of the initial topography as a function
         * of surface position.
         */
        double
        value (const Point<dim-1> &p) const override;

        /**
         * Return the maximum value of the elevation.
         */
        double max_topography () const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * The maximum value the topography can take.
         *
         * This variable is read from the parameter file through a parameter called 'Maximum topography value'.
         */
        double max_topo;

        /**
         * A function object representing the topography.
         */
        Functions::ParsedFunction<dim> initial_topography_function;

        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are cartesian and spherical.
         *
         * This variable is read from the parameter file through a parameter called 'Coordinate system'.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;
    };
  }
}


#endif
