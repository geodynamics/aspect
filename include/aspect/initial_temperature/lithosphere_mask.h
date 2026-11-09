// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2019-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_initial_temperature_lithosphere_mask_h
#define _aspect_initial_temperature_lithosphere_mask_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    namespace LABDepth
    {
      template <int dim>
      class LABDepthLookup : public SimulatorAccess<dim>
      {
        public:
          /**
           * Empty Constructor.
           */
          LABDepthLookup ();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program. Checks preconditions.
           */
          void
          initialize ();

          /**
           * Return LAB depth as a function of position (latitude and longitude). This
           * function returns either a constant value or values from a text file,
           * depending on the input parameters. Text files are read in
           * two dimensions so the third column (depth) is treated as data.
           */
          double
          get_lab_depth (const Point<dim> &position) const;

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
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * Reads in file containing input data.
           */
          Utilities::StructuredDataLookup<2> lab_depths;

          /**
           * Directory in which the LAB depth file is present.
           *
           * This variable is read from the parameter file through a parameter called 'Data directory'.
           */
          std::string data_directory;

          /**
           * File name of the LAB depth file.
           * This variable is read from the parameter file through a parameter called 'LAB depth filename'.
           */
          std::string LAB_file_name;

          /**
           * An enum to describe where the LAB depth is coming from.
           */
          enum LABDepthSource
          {
            Value,
            File
          };

          /**
           * Currently chosen source for the LAB depth.
           */
          LABDepthSource LAB_depth_source;

          /**
           * This parameter gives the maximum depth of the lithosphere. The
           * model returns nans below this depth. This parameter is only used if LAB
           * depth source is set to 'Value'.
           * This variable is read from the parameter file through a parameter called 'Maximum lithosphere depth'.
           */
          double max_depth;
      };
    }

    /**
     * A class that implements a constant reference temperature
     * above a specified depth and nans below the specified depth.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class LithosphereMask : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * Return the initial temperature as a function of position. For the
         * current class, this function uses the given lithosphere temperature
         * above the lithosphere-asthenosphere boundary and nans below.
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
         * This parameter gives the initial temperature set within the lithosphere.
         *
         * This variable is read from the parameter file through a parameter called 'Lithosphere temperature'.
         */
        double lithosphere_temperature;

        /**
         *
         */
        LABDepth::LABDepthLookup<dim> lab_depth_lookup;
    };

  }
}

#endif
