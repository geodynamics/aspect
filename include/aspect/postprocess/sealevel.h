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


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs the sea level to file.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class SeaLevel : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Sets up structured data lookup for topography and ice height input data.
         */
        void initialize() override;

        /**
         * Output non-uniform sea level change [m] to file.
         */
        std::pair<std::string,std::string> 
        execute (TableHandler &statistics) override;

        /**
         * Evaluate the total surface pressure from ice and water loading at 
         * a point.
         */
        double
        compute_total_surface_pressure(const Point<dim> &) const;

        /**
         * Register with the simulator the other postprocessors that we need
         * (namely: geoid).
         */
        std::list<std::string>
        required_other_postprocessors() const override;

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
         * Function to determine non-uniform sea level change for any given location.
         */
        double
        compute_nonuniform_sealevel_change(const Point<dim> &position) const;
        
        /**
         * Function to determine the sea level offset.
         */
        double
        compute_sealevel_offset(const double &outer_radius) const;

        /**
         * Information about the location of topography data files.
         */
        std::string data_directory_topography;
        std::string data_file_name_topography;

        /**
         * Information about the location of ice height data files.
         */        
        std::string data_directory_iceheight;
        std::string data_file_name_iceheight;

        /**
         * Pointer to the StructuredDataLookup object that holds the topography  and ice height data.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<2>> topography_lookup;
        std::unique_ptr<Utilities::StructuredDataLookup<2>> iceheight_lookup;

        /**
         * description
         */
        double density_water;

        /**
         * description
         */
        double density_ice;

        /**
         * The sealevel offset.
         */
        double sealevel_offset;

        /**
         * Whether or not to produce text files with sea level values.
         */
        bool write_to_file;

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
