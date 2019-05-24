/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_lithosphere_mask_h
#define _aspect_initial_temperature_lithosphere_mask_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

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
         * Empty Constructor.
         */
        LithosphereMask ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        /**
         * Directory in which the LAB depth file is present.
         */
        std::string data_directory;

        /**
         * File name of the LAB depth file.
         */
        std::string LAB_file_name;

        /**
         * Return LAB depth as a function of position (latitude and longitude). For the
         * current class, this function returns value from the text files.
         */
        double
        ascii_lab (const Point<2> &position) const;

        /**
         * This parameter gives the maximum depth of the lithosphere. The
         * model be nans below this depth.
         */
        double max_depth;

        /**
         * This parameter gives the initial temperature set within the lithosphere.
         */
        double lithosphere_temperature;

        /**
         * Return the initial temperature as a function of position. For the
         * current class, this function uses the given lithosphere temperature
         * above max_depth and nans below max_depth.
         */
        double
        initial_temperature (const Point<dim> &position) const;

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
        Utilities::AsciiDataLookup<2> lab_depths;

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
    };
  }
}


#endif
