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

#ifndef _aspect_postprocess_particle_h
#define _aspect_postprocess_particle_h

#include <aspect/postprocess/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/particle/property/interface.h>

#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/data_out_base.h>
#include <tuple>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      /**
       * This class is responsible for writing the particle data into a format that can
       * be written by deal.II, in particular a list of 'patches' that contain one
       * particle per patch. The base class DataOutInterface is templated with a
       * dimension of zero (the dimension of the particle / point), and a space dimension
       * of dim (the dimension in which this zero-dimensional particle lives).
       */
      template<int dim>
      class ParticleOutput : public dealii::DataOutInterface<0,dim>
      {
        public:
          /**
           * This function prepares the data for writing. It reads the data from @p particle_handler and their
           * property information from @p property_information, and builds a list of patches that is stored
           * internally until the destructor is called. This function needs to be called before one of the
           * write function of the base class can be called to write the output data.
           */
          void build_patches(const Particles::ParticleHandler<dim> &particle_handler,
                             const aspect::Particle::Property::ParticlePropertyInformation &property_information,
                             const std::vector<std::string> &exclude_output_properties,
                             const bool only_group_3d_vectors);

        private:
          /**
           * Implementation of the corresponding function of the base class.
           */
          const std::vector<DataOutBase::Patch<0,dim> > &
          get_patches () const override;

          /**
           * Implementation of the corresponding function of the base class.
           */
          std::vector< std::string >
          get_dataset_names () const override;

          /**
           * Implementation of the corresponding function of the base class.
           */
#if DEAL_II_VERSION_GTE(9,1,0)
          std::vector<
          std::tuple<unsigned int,
              unsigned int,
              std::string,
              DataComponentInterpretation::DataComponentInterpretation> >
              get_nonscalar_data_ranges () const override;
#else
          std::vector<std::tuple<unsigned int, unsigned int, std::string> >
          get_vector_data_ranges() const override;
#endif

          /**
           * Output information that is filled by build_patches() and
           * written by the write function of the base class.
           */
          std::vector<DataOutBase::Patch<0,dim> > patches;

          /**
           * A list of field names for all data components stored in patches.
           */
          std::vector<std::string> dataset_names;

          /**
           * Store which of the data fields are vectors.
           */
#if DEAL_II_VERSION_GTE(9,1,0)
          std::vector<
          std::tuple<unsigned int,
              unsigned int,
              std::string,
              DataComponentInterpretation::DataComponentInterpretation> >
              vector_datasets;
#else
          std::vector<std::tuple<unsigned int, unsigned int, std::string> >
          vector_datasets;
#endif
      };
    }

    /**
     * A Postprocessor that creates particles, which follow the
     * velocity field of the simulation. The particles can be generated
     * and propagated in various ways and they can carry a number of
     * constant or time-varying properties. The postprocessor can write
     * output positions and properties of all particles at chosen intervals,
     * although this is not mandatory. It also allows other parts of the
     * code to query the particles for information.
     */
    template <int dim>
    class Particles : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Particles();

        /**
         * Destructor.
         */
        ~Particles() override;

        /**
         * Execute this postprocessor. Derived classes will implement this
         * function to do whatever they want to do to evaluate the solution at
         * the current time step.
         *
         * @param[in,out] statistics An object that contains statistics that
         * are collected throughout the simulation and that will be written to
         * an output file at the end of each time step. Postprocessors may
         * deposit data in these tables for later visualization or further
         * processing.
         *
         * @return A pair of strings that will be printed to the screen after
         * running the postprocessor in two columns; typically the first
         * column contains a description of what the data is and the second
         * contains a numerical value of this data. If there is nothing to
         * print, simply return two empty strings.
         */
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

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
         * to create a master file that can make the association
         * between simulation time and corresponding file name (this
         * is done because there is no way to store the simulation
         * time inside the .pvtu or .vtu files).
         */
        std::vector<std::pair<double,std::string> > times_and_pvtu_file_names;

        /**
         * A corresponding variable that we use for the .visit files created
         * by DataOutInterface::write_visit_record. The second part of a
         * pair contains all files that together form a time step.
         */
        std::vector<std::pair<double,std::vector<std::string> > > times_and_vtu_file_names;

        /**
         * A list of list of filenames, sorted by timestep, that correspond to
         * what has been created as output. This is used to create a master
         * .visit file for the entire simulation.
         */
        std::vector<std::vector<std::string> > output_file_names_by_timestep;

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
         * Handle to a thread that is used to write data in the background.
         * The writer() function runs on this background thread.
         */
        Threads::Thread<void> background_thread;

        /**
         * Stores the particle property fields which are excluded from output
         * to the visualization file.
         */
        std::vector<std::string> exclude_output_properties;

        /**
         * A function that writes the text in the second argument to a file
         * with the name given in the first argument. The function is run on a
         * separate thread to allow computations to continue even though
         * writing data is still continuing. The function takes over ownership
         * of these arguments and deletes them at the end of its work.
         */
        static
        void writer (const std::string filename,
                     const std::string temporary_filename,
                     const std::string *file_contents);

        /**
         * Write the various master record files. The master files are used by
         * visualization programs to identify which of the output files in a
         * directory, possibly one file written by each processor, belong to a
         * single time step and/or form the different time steps of a
         * simulation. For Paraview, this is a <code>.pvtu</code> file per
         * time step and a <code>.pvd</code> for all time steps. For Visit it
         * is a <code>.visit</code> file per time step and one for all time
         * steps.
         *
         * @param data_out The DataOut object that was used to write the
         * solutions.
         * @param solution_file_prefix The stem of the filename to be written.
         * @param filenames List of filenames for the current output from all
         * processors.
         */
        void write_master_files (const internal::ParticleOutput<dim> &data_out,
                                 const std::string &solution_file_prefix,
                                 const std::vector<std::string> &filenames);
    };
  }
}

#endif
