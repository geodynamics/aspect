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


#ifndef _aspect_initial_temperature_patch_on_S40RTS_h
#define _aspect_initial_temperature_patch_on_S40RTS_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_temperature/S40RTS_perturbation.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that implements a prescribed temperature field determined from
     * an upper mantle Vs model (input as an ascii file) above a specified depth
     * and S40RTS below the specified depth.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class PatchOnS40RTS : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        PatchOnS40RTS ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * Return Vs as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        ascii_grid_vs (const Point<dim> &position) const;

        /**
         * This parameter gives the maximum depth of the Vs ascii grid. The
         * model will read in Vs from S40RTS below this depth.
         *
         * This variable is read from the parameter file through a parameter called 'Maximum grid depth'.
         */
        double max_grid_depth;

        /**
         * This parameter gives the range (above maximum grid depth) over which to smooth.
         * Smoothing is done with a depth weighted combination of the values in the ascii grid and S40RTS
         * at each point.
         * This variable is read from the parameter file through a parameter called 'Smoothing length scale'.
         */
        double smoothing_length_scale;

        /**
         * Return the initial temperature as a function of position. For the
         * current class, this function calculates temperature from ascii grid Vs data
         * above max_grid_depth and S40RTS Vs data below max_grid_depth.
         */
        double
        initial_temperature (const Point<dim> &position) const override;

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
         * Declare a member variable of type S40RTSPerturbation that allows us to call
         * functions from S40RTS_perturbation.cc.
         */
        S40RTSPerturbation<dim> s40rts;

        /**
         * This parameter is the depth down to which shear wave perturbations are
         * zeroed out.
         *
         * This variable is read from the parameter file through a parameter called 'Remove temperature heterogeneity down to specified depth'.
         */
        double no_perturbation_depth_patch;

    };
  }
}


#endif
