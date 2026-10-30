// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_adiabatic_conditions_ascii_data_h
#define _aspect_adiabatic_conditions_ascii_data_h


#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/point.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    /**
     * A simple class that reads adiabatic conditions from a file.
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataProfile<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        AsciiData ();

        /**
         * Initialization function. Because this function is called after
         * initializing the SimulatorAccess, all of the necessary information
         * is available to calculate the adiabatic profile. It computes the
         * adiabatic conditions along a vertical transect of the geometry
         * based on the given material model and other quantities.
         */
        void initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataProfile<dim>::initialize;

        /**
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated. This is for example the case for the simple
         * compressible material model, which uses the adiabatic temperature
         * as reference temperature to calculate the density. For the
         * calculation of the adiabatic conditions this functionality is
         * simply switched off, because we are always on the reference
         * profile. This way the plugin behaves differently at initialization
         * time of the adiabatic conditions and during the main model run.
         */
        bool is_initialized() const override;

        /**
         * Return the adiabatic temperature at a given point of the domain.
         */
        double temperature (const Point<dim> &p) const override;

        /**
         * Return the adiabatic pressure at a given point of the domain.
         */
        double pressure (const Point<dim> &p) const override;

        /**
         * Return the reference density at a given point of the domain.
         */
        double density (const Point<dim> &p) const override;

        /**
         * Return the derivative of the density with respect to depth
         * at the given point @p p.
         */
        double density_derivative (const Point<dim> &p) const override;

        /**
         * Declare the parameters for the input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * Whether the adiabatic conditions are already calculated. This is
         * important for plugins that are used by the adiabatic conditions but
         * also depend on the adiabatic conditions. This way they can behave
         * differently in initialization and model run.
         */
        bool initialized;

        /**
         * The column indices of the temperature, pressure, and density column
         * in the data file.
         */
        unsigned int temperature_index;
        unsigned int pressure_index;
        unsigned int density_index;
    };
  }
}


#endif
