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


#ifndef _aspect_postprocess_sealevel_h
#define _aspect_postprocess_sealevel_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/initial_temperature/interface.h>

namespace aspect
{

  
  namespace Postprocess
  {
    template <int dim>
    class SeaLevel : public Interface<dim>, public ::aspect::SimulatorAccess<dim> //, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Output sea level [m] to file
         */
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        //

        /**
         * Register with the simulator the other postprocessors that we need
         * (namely: topogaphy, geoid, selfgravitation).
         */
        std::list<std::string>
        required_other_postprocessors() const override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
        /**
         * @}
         */

      private:
        /**
         * Function to determine sea level for any given location
         */
        double
        sea_level_equation(const Point<dim> &position);

        /**
         * Whether or not to produce text files with sea level values
         */
        bool write_to_file;

        /**
         * Pointer to the structured data
         */
        std::unique_ptr<Utilities::StructuredDataLookup<dim-1>> data_lookup;
        /**
         * Directory in which the input data files, i.e., GPlates model or ascii data
         * are present.
         */
        std::string data_directory;

        /**
         * Filename of the input Gplates model or ascii data file. For GPlates, the file names
         * can contain the specifiers %s and/or %c (in this order), meaning the name of the
         * boundary and the number of the data file time step.
         */
        std::string data_file_name;

        /**
         * Interval between the generation of text output. This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last text output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        double last_output_time;

    };
  }
}


#endif
