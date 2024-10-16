/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_sea_level_h
#define _aspect_postprocess_sea_level_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs the sea level to file. The sea level is based
     * on the ice heights input data, the solid Earth topography input data for
     * defining the ocean basin, the geoid displacement from the geoid postprocessor,
     * and the free surface topography from the geometry model. The sea level
     * postprocessor requires a 3D spherical shell geometry with a free surface.
     * The sea level computation is based on \\cite{Martinec2018}.
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
         * a point at the top surface of the domain. While the function will
         * return pressure values for any location, it will always assume the
         * volume between the point and the sea level is filled with water or ice,
         * therefore calling this function for positions inside the domain will result
         * in unrealistic pressure values.
         */
        double
        compute_total_surface_pressure(const Point<dim> &position) const;

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

        /**
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * Function to determine non-uniform sea level change for any given location.
         * This is a local computation that makes use of the computed sea level offset,
         * so make sure to call compute_sea_level_offset() before calling this function.
         */
        double
        compute_nonuniform_sea_level_change(const Point<dim> &position) const;

        /**
         * Function to determine the global sea level offset.
         */
        double
        compute_sea_level_offset();

        /**
         * Information about the location of topography data files.
         */
        std::string data_directory_topography;
        std::string data_file_name_topography;

        /**
         * Information about the location of ice height data files.
         */
        std::string data_directory_ice_height;
        std::string data_file_name_ice_height;

        /**
         * Pointer to the StructuredDataLookup object that holds the topography and ice height data.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<2>> topography_lookup;
        std::unique_ptr<Utilities::StructuredDataLookup<2>> ice_height_lookup;

        /**
         * The density of water.
         */
        double density_water;

        /**
         * The density of ice.
         */
        double density_ice;

        /**
         * The global sea level offset.
         */
        double sea_level_offset;

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
