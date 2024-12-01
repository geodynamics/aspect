/*
 Copyright (C) 2023 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_postprocess_crystal_preferred_orientation_h
#define _aspect_postprocess_crystal_preferred_orientation_h

#include <aspect/postprocess/interface.h>
#include <aspect/particle/manager.h>

#include <aspect/simulator_access.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/data_out_base.h>
#include <tuple>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A Postprocessor that writes out CPO specific particle data.
     * It can write out the CPO data as it is stored (raw) and/or as a
     * random draw volume weighted representation. The latter one
     * is recommended for plotting against real data. For both representations
     * the specific output fields and their order can be set.
     */
    template <int dim>
    class CrystalPreferredOrientation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        CrystalPreferredOrientation();

        /**
         * Destructor.
         */
        ~CrystalPreferredOrientation();


        /**
         * Initialize function.
         */
        void initialize () override;

        /**
         * A Postprocessor that writes out CPO specific particle data.
         * It can write out the CPO data as it is stored (raw) and/or as a
         * random draw volume weighted representation. The latter one
         * is recommended for plotting against real data. For both representations
         * the specific output fields and their order can be set.
         */
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        /**
         * This function ensures that the particle postprocessor is run before
         * this postprocessor.
         */
        std::list<std::string>
        required_other_postprocessors () const override;

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
         * Stores the simulation end time, so that it always produces output
         * at the last timestep.
         */
        double end_time;

        /**
         * Enums specifying what information to write:
         *
         * VolumeFraction: Write the volume fraction
         * RotationMatrix: Write the whole rotation matrix of a grain (in deal.ii order)
         * EulerAngles: Write out the Z-X-Z Euler angle of a grain
         * not_found: Error enum, the requested output was not found
         */
        enum class Output
        {
          VolumeFraction, RotationMatrix, EulerAngles, not_found
        };

        /**
         * Converts a string to an Output enum.
         */
        Output string_to_output_enum(const std::string &string);

        /**
         * Random number generator used for random draw volume weighting.
         */
        mutable std::mt19937 random_number_generator;

        /**
         * Random number generator seed used to initialize the random number generator.
         */
        unsigned int random_number_seed;

        /**
         * Interval between output (in years if appropriate simulation
         * parameter is set, otherwise seconds)
         */
        double output_interval;

        /**
         * Records time for next output to occur
         */
        double last_output_time;

        /**
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);

        /**
         * Consecutively counted number indicating the how-manyth time we will
         * create output the next time we get to it.
         */
        unsigned int output_file_number;

        /**
         * Graphical output format.
         */
        std::vector<std::string> output_formats;

        /**
         * A list of pairs (time, pvtu_filename) that have so far been written
         * and that we will pass to DataOutInterface::write_pvd_record
         * to create a main file that can make the association
         * between simulation time and corresponding file name (this
         * is done because there is no way to store the simulation
         * time inside the .pvtu or .vtu files).
         */
        std::vector<std::pair<double,std::string>> times_and_pvtu_file_names;

        /**
         * A corresponding variable that we use for the .visit files created
         * by DataOutInterface::write_visit_record. The second part of a
         * pair contains all files that together form a time step.
         */
        std::vector<std::pair<double,std::vector<std::string>>> times_and_vtu_file_names;

        /**
         * A list of list of filenames, sorted by timestep, that correspond to
         * what has been created as output. This is used to create a main
         * .visit file for the entire simulation.
         */
        std::vector<std::vector<std::string>> output_file_names_by_timestep;

        /**
         * A set of data related to XDMF file sections describing the HDF5
         * heavy data files created. These contain things such as the
         * dimensions and names of data written at all steps during the
         * simulation.
         */
        std::vector<XDMFEntry>  xdmf_entries;

        /**
         * VTU file output supports grouping files from several CPUs into one
         * file using MPI I/O when writing on a parallel filesystem. 0 means
         * no grouping (and no parallel I/O). 1 will generate one big file
         * containing the whole solution.
         */
        unsigned int group_files;

        /**
         * On large clusters it can be advantageous to first write the
         * output to a temporary file on a local file system and later
         * move this file to a network file system. If this variable is
         * set to a non-empty string it will be interpreted as a temporary
         * storage location.
         */
        std::string temporary_output_location;

        /**
         * File operations can potentially take a long time, blocking the
         * progress of the rest of the model run. Setting this variable to
         * 'true' moves this process into a background thread, while the
         * rest of the model continues.
         */
        bool write_in_background_thread;

        /**
         * Handle to a thread that is used to write main file data in the
         * background. The writer() function runs on this background thread.
         */
        std::thread background_thread_main;

        /**
         * What "raw" CPO data to write out.
         */
        std::vector<std::pair<unsigned int,Output>> write_raw_cpo;

        /**
         * Whether computing raw Euler angles is needed.
         */
        bool compute_raw_euler_angles;

        /**
         * Handle to a thread that is used to write content file data in the
         * background. The writer() function runs on this background thread.
         */
        std::thread background_thread_content_raw;

        /**
         * What "draw volume weighted" CPO data to write out. Draw volume weighted means
         * the grain properties (size and rotation) are put in a list sorted based by volume
         * and picked randomly, with large volumes having a higher chance of being picked.
         */
        std::vector<std::pair<unsigned int,Output>> write_draw_volume_weighted_cpo;

        /**
         * Whether computing weighted A matrix is needed.
         */
        bool compute_weighted_rotation_matrix;

        /**
         * Handle to a thread that is used to write content file data in the
         * background. The writer() function runs on this background thread.
         */
        std::thread background_thread_content_draw_volume_weighting;

        /**
         * Whether to compress the raw and weighed cpo data output files with zlib.
         */
        bool compress_cpo_data_files;

        /**
         * A function that writes the text in the second argument to a file
         * with the name given in the first argument. The function is run on a
         * separate thread to allow computations to continue even though
         * writing data is still continuing. The function takes over ownership
         * of these arguments and deletes them at the end of its work.
         */
        static
        void writer (const std::string &filename,
                     const std::string &temporary_filename,
                     const std::string &file_contents,
                     const bool compress_contents);
    };
  }
}

#endif
