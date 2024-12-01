/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/manager.h>
#include <aspect/utilities.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cstdio>
#include <unistd.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      template <int dim>
      void
      ParticleOutput<dim>::build_patches(const dealii::Particles::ParticleHandler<dim> &particle_handler,
                                         const aspect::Particle::Property::ParticlePropertyInformation &property_information,
                                         const std::vector<std::string> &exclude_output_properties,
                                         const bool only_group_3d_vectors)
      {
        // First store the names of the data fields that should be written
        dataset_names.reserve(property_information.n_components()+1);
        dataset_names.emplace_back("id");

        // This is a map from an index for the particle property vector to an index in the output
        // vector. Values equal 0 indicate a property will not be written, values bigger 0
        // indicate into which output entry the property value should be written.
        std::vector<unsigned int> property_index_to_output_index(property_information.n_components(),0);

        if (std::find(exclude_output_properties.begin(),
                      exclude_output_properties.end(),
                      "all")
            == exclude_output_properties.end())
          {
            // Start the output index from 1, because 0 is occupied by the "id"
            unsigned int output_index = 1;
            for (unsigned int field_index = 0; field_index < property_information.n_fields(); ++field_index)
              {
                // Determine some info about the field.
                const unsigned n_components = property_information.get_components_by_field_index(field_index);
                const std::string field_name = property_information.get_field_name_by_index(field_index);
                const unsigned int field_position = property_information.n_fields() == 0
                                                    ?
                                                    0
                                                    :
                                                    property_information.get_position_by_field_index(field_index);

                // HDF5 only supports 3d vector output, therefore only treat output fields as vector if we
                // have a dimension of 3 and 3 components.
                const bool field_is_vector = (only_group_3d_vectors == false
                                              ?
                                              n_components == dim
                                              :
                                              dim == 3 && n_components == 3);

                // Determine if this field should be excluded, if so, skip it
                bool exclude_property = false;
                for (const auto &property : exclude_output_properties)
                  if (field_name.find(property) != std::string::npos)
                    {
                      exclude_property = true;
                      break;
                    }

                if (exclude_property == false)
                  {
                    // For each component record its name and position in output vector
                    for (unsigned int component_index=0; component_index<n_components; ++component_index)
                      {
                        // If it is a 1D element, or a vector, print just the name, otherwise append the index after an underscore
                        if ((n_components == 1) || field_is_vector)
                          dataset_names.push_back(field_name);
                        else
                          dataset_names.push_back(field_name + "_" + Utilities::to_string(component_index));

                        property_index_to_output_index[field_position + component_index] = output_index;
                        ++output_index;
                      }

                    // If the property has dim components, we treat it as vector
                    if (n_components == dim)
                      {
                        vector_datasets.emplace_back(property_index_to_output_index[field_position],
                                                     property_index_to_output_index[field_position]+n_components-1,
                                                     field_name,
                                                     DataComponentInterpretation::component_is_part_of_vector);
                      }
                  }
              }
          }

        // Now build the actual patch data
        patches.resize(particle_handler.n_locally_owned_particles());
        typename dealii::Particles::ParticleHandler<dim>::particle_iterator particle = particle_handler.begin();

        for (types::particle_index i=0; particle != particle_handler.end(); ++particle, ++i)
          {
            patches[i].vertices[0] = particle->get_location();
            patches[i].patch_index = i;

            patches[i].data.reinit(dataset_names.size(),1);

            patches[i].data(0,0) = particle->get_id();

            if (particle->has_properties())
              {
                const ArrayView<const double> properties = particle->get_properties();

                for (unsigned int property_index = 0; property_index < properties.size(); ++property_index)
                  {
                    if (property_index_to_output_index[property_index] > 0)
                      patches[i].data(property_index_to_output_index[property_index],0) = properties[property_index];
                  }
              }
          }
      }

      template <int dim>
      const std::vector<DataOutBase::Patch<0,dim>> &
      ParticleOutput<dim>::get_patches () const
      {
        return patches;
      }

      template <int dim>
      std::vector<std::string>
      ParticleOutput<dim>::get_dataset_names () const
      {
        return dataset_names;
      }

      template <int dim>
      std::vector<
      std::tuple<unsigned int,
          unsigned int, std::string,
          DataComponentInterpretation::DataComponentInterpretation>>
          ParticleOutput<dim>::get_nonscalar_data_ranges () const
      {
        return vector_datasets;
      }
    }

    template <int dim>
    Particles<dim>::Particles ()
      :
      // the following value is later read from the input file
      output_interval (0),
      // initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      last_output_time (std::numeric_limits<double>::quiet_NaN())
      ,output_file_number (numbers::invalid_unsigned_int),
      group_files(0),
      write_in_background_thread(false)
    {}



    template <int dim>
    Particles<dim>::~Particles ()
    {
      // make sure a thread that may still be running in the background,
      // writing data, finishes
      if (background_thread.joinable())
        background_thread.join ();
    }



    template <int dim>
    // We need to pass the arguments by value, as this function can be called on a separate thread:
    void Particles<dim>::writer (const std::string &filename, //NOLINT(performance-unnecessary-value-param)
                                 const std::string &temporary_output_location, //NOLINT(performance-unnecessary-value-param)
                                 const std::string &file_contents)
    {
      std::string tmp_filename = filename;
      if (temporary_output_location != "")
        {
          tmp_filename = temporary_output_location + "/aspect.tmp.XXXXXX";

          // Create the temporary file; get at the actual filename
          // by using a C-style string that mkstemp will then overwrite
          std::vector<char> tmp_filename_x (tmp_filename.size()+1);
          std::strcpy(tmp_filename_x.data(), tmp_filename.c_str());
          const int tmp_file_desc = mkstemp(tmp_filename_x.data());
          tmp_filename = tmp_filename_x.data();

          // If we failed to create the temp file, just write directly to the target file.
          // We also provide a warning about this fact. There are places where
          // this fails *on every node*, so we will get a lot of warning messages
          // into the output; in these cases, just writing multiple pieces to
          // std::cerr will produce an unreadable mass of text; rather, first
          // assemble the error message completely, and then output it atomically
          if (tmp_file_desc == -1)
            {
              const std::string x = ("***** WARNING: could not create temporary file <"
                                     +
                                     tmp_filename
                                     +
                                     ">, will output directly to final location. This may negatively "
                                     "affect performance. (On processor "
                                     + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD))
                                     + ".)\n");

              std::cerr << x << std::flush;

              tmp_filename = filename;
            }
          else
            close(tmp_file_desc);
        }

      std::ofstream out(tmp_filename);

      AssertThrow (out, ExcMessage(std::string("Trying to write to file <") +
                                   filename +
                                   ">, but the file can't be opened!"));

      // now write and then move the tmp file to its final destination
      // if necessary
      out << file_contents;
      out.close ();

      if (tmp_filename != filename)
        {
          std::string command = std::string("mv ") + tmp_filename + " " + filename;
          int error = system(command.c_str());

          AssertThrow(error == 0,
                      ExcMessage("Could not move " + tmp_filename + " to "
                                 + filename + ". On processor "
                                 + Utilities::int_to_string(Utilities::MPI::this_mpi_process (MPI_COMM_WORLD)) + "."));
        }
    }

    template <int dim>
    void
    Particles<dim>::write_description_files (const internal::ParticleOutput<dim> &data_out,
                                             const std::string &description_file_prefix,
                                             const std::string &solution_file_directory,
                                             const std::string &solution_file_prefix,
                                             const std::vector<std::string> &filenames)
    {
      const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                               this->get_time() / year_in_seconds :
                                               this->get_time());
      const std::string pvtu_filename = (solution_file_prefix +
                                         ".pvtu");
      std::ofstream pvtu_file (this->get_output_directory() + solution_file_directory +
                               pvtu_filename);
      data_out.write_pvtu_record (pvtu_file, filenames);

      // now also generate a .pvd file that matches simulation
      // time and corresponding .pvtu record
      times_and_pvtu_file_names[description_file_prefix].emplace_back(time_in_years_or_seconds, solution_file_directory + pvtu_filename);

      const std::string pvd_filename = (this->get_output_directory() + description_file_prefix + ".pvd");
      std::ofstream pvd_file (pvd_filename);

      DataOutBase::write_pvd_record (pvd_file, times_and_pvtu_file_names[description_file_prefix]);

      // finally, do the same for VisIt via the .visit file for this
      // time step, as well as for all time steps together
      const std::string visit_filename = (this->get_output_directory()
                                          + solution_file_directory
                                          + solution_file_prefix
                                          + ".visit");
      std::ofstream visit_file (visit_filename);

      DataOutBase::write_visit_record (visit_file, filenames);

      {
        // the global .visit file needs the relative path because it sits a
        // directory above
        std::vector<std::string> filenames_with_path;
        filenames_with_path.reserve(filenames.size());
        for (const auto &filename : filenames)
          {
            filenames_with_path.push_back(solution_file_directory + filename);
          }

        output_file_names_by_timestep[description_file_prefix].push_back (filenames_with_path);
      }

      std::ofstream global_visit_file (this->get_output_directory() +
                                       description_file_prefix + ".visit");

      std::vector<std::pair<double, std::vector<std::string>>> times_and_output_file_names;
      for (unsigned int timestep=0; timestep<times_and_pvtu_file_names[description_file_prefix].size(); ++timestep)
        times_and_output_file_names.emplace_back(times_and_pvtu_file_names[description_file_prefix][timestep].first,
                                                 output_file_names_by_timestep[description_file_prefix][timestep]);
      DataOutBase::write_visit_record (global_visit_file, times_and_output_file_names);
    }

    template <int dim>
    std::pair<std::string,std::string>
    Particles<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        last_output_time = this->get_time() - output_interval;

      bool write_output = true;
      std::string number_of_advected_particles = "";
      std::string screen_output = "";

      for (unsigned int particle_manager = 0; particle_manager < this->n_particle_managers(); ++particle_manager)
        {
          const Particle::Manager<dim> &manager = this->get_particle_manager(particle_manager);

          const std::string statistics_column_name = (particle_manager == 0 ?
                                                      "Number of advected particles" :
                                                      "Number of advected particles (Particle system " + Utilities::int_to_string(particle_manager+1) + ")");

          statistics.add_value(statistics_column_name,manager.n_global_particles());

          if (particle_manager > 0)
            {
              number_of_advected_particles += ", ";
            }
          number_of_advected_particles += Utilities::int_to_string(manager.n_global_particles());
        }

      // If it's not time to generate an output file
      // return early with the number of particles that were advected
      if (this->get_time() < last_output_time + output_interval)
        {
          write_output = false;
        }

      // If we do not write output
      // return early with the number of particles that were advected
      if (output_formats.size() == 0 || output_formats[0] == "none")
        {
          // Up the next time we need output. This is relevant to correctly
          // write output after a restart if the format is changed.
          set_last_output_time (this->get_time());

          write_output = false;
        }

      if (write_output == false)
        return std::make_pair("Number of advected particles:", number_of_advected_particles);

      if (output_file_number == numbers::invalid_unsigned_int)
        output_file_number = 0;
      else
        ++output_file_number;

      for (unsigned int particle_manager = 0; particle_manager < this->n_particle_managers(); ++particle_manager)
        {
          const Particle::Manager<dim> &manager = this->get_particle_manager(particle_manager);
          std::string particles_output_base_name = "particles";
          if (particle_manager > 0)
            {
              particles_output_base_name += "-" + Utilities::int_to_string(particle_manager+1);
            }

          // Create the particle output
          const bool output_hdf5 = std::find(output_formats.begin(), output_formats.end(),"hdf5") != output_formats.end();
          internal::ParticleOutput<dim> data_out;
          data_out.build_patches(manager.get_particle_handler(),
                                 manager.get_property_manager().get_data_info(),
                                 exclude_output_properties,
                                 output_hdf5);

          // Now prepare everything for writing the output and choose output format
          std::string particle_file_prefix = particles_output_base_name + "-" + Utilities::int_to_string (output_file_number, 5);

          const double time_in_years_or_seconds = (this->convert_output_to_years() ?
                                                   this->get_time() / year_in_seconds :
                                                   this->get_time());

          for (const auto &output_format : output_formats)
            {
              // this case was handled above
              Assert(output_format != "none", ExcInternalError());

              if (output_format == "hdf5")
                {
                  const std::string particle_file_name = particles_output_base_name + "/" + particle_file_prefix + ".h5";
                  const std::string xdmf_filename = particles_output_base_name + ".xdmf";

                  // Do not filter redundant values, there are no duplicate particles
                  DataOutBase::DataOutFilter data_filter(DataOutBase::DataOutFilterFlags(false, true));

                  data_out.write_filtered_data(data_filter);
                  data_out.write_hdf5_parallel(data_filter,
                                               this->get_output_directory()+particle_file_name,
                                               this->get_mpi_communicator());

                  const XDMFEntry new_xdmf_entry = data_out.create_xdmf_entry(data_filter,
                                                                              particle_file_name,
                                                                              time_in_years_or_seconds,
                                                                              this->get_mpi_communicator());
                  xdmf_entries[particles_output_base_name].push_back(new_xdmf_entry);
                  data_out.write_xdmf_file(xdmf_entries[particles_output_base_name], this->get_output_directory() + xdmf_filename,
                                           this->get_mpi_communicator());
                }
              else if (output_format == "vtu")
                {
                  // Write descriptive files (.pvtu,.pvd,.visit) on the root process
                  const int my_id = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

                  if (my_id == 0)
                    {
                      std::vector<std::string> filenames;
                      const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());
                      const unsigned int n_files = (group_files == 0) ? n_processes : std::min(group_files,n_processes);
                      for (unsigned int i=0; i<n_files; ++i)
                        filenames.push_back (particle_file_prefix
                                             + "." + Utilities::int_to_string(i, 4)
                                             + ".vtu");
                      write_description_files (data_out,
                                               particles_output_base_name,
                                               particles_output_base_name + "/",
                                               particle_file_prefix,
                                               filenames);
                    }

                  const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

                  const unsigned int my_file_id = (group_files == 0
                                                   ?
                                                   my_id
                                                   :
                                                   my_id % group_files);
                  const std::string filename = this->get_output_directory()
                                               + particles_output_base_name
                                               + "/"
                                               + particle_file_prefix
                                               + "."
                                               + Utilities::int_to_string (my_file_id, 4)
                                               + ".vtu";

                  // pass time step number and time as metadata into the output file
                  DataOutBase::VtkFlags vtk_flags;
                  vtk_flags.cycle = this->get_timestep_number();
                  vtk_flags.time = time_in_years_or_seconds;

                  data_out.set_flags (vtk_flags);

                  // Write as many files as processes. For this case we support writing in a
                  // background thread and to a temporary location, so we first write everything
                  // into a string that is written to disk in a writer function
                  if ((group_files == 0) || (group_files >= n_processes))
                    {
                      // Put the content we want to write into a string object that
                      // we can then write in the background
                      std::unique_ptr<std::string> file_contents;
                      {
                        std::ostringstream tmp;

                        data_out.write (tmp, DataOutBase::parse_output_format(output_format));
                        file_contents = std::make_unique<std::string>(tmp.str());
                      }

                      if (write_in_background_thread)
                        {
                          // Wait for all previous write operations to finish, should
                          // any be still active, ...
                          if (background_thread.joinable())
                            background_thread.join ();

                          // ...then continue with writing our own data.
                          background_thread
                            = std::thread([ my_filename = filename, // filename is const, so we can not move from it
                                            my_temporary_output_location = temporary_output_location,
                                            my_file_contents = std::move(file_contents)]()
                          {
                            writer (my_filename, my_temporary_output_location, *my_file_contents);
                          });
                        }
                      else
                        writer(filename,temporary_output_location,*file_contents);
                    }
                  // Just write one data file in parallel
                  else if (group_files == 1)
                    {
                      data_out.write_vtu_in_parallel(filename,
                                                     this->get_mpi_communicator());
                    }
                  // Write as many output files as 'group_files' groups
                  else
                    {
                      int color = my_id % group_files;

                      MPI_Comm comm;
                      int ierr = MPI_Comm_split(this->get_mpi_communicator(), color, my_id, &comm);
                      AssertThrowMPI(ierr);

                      data_out.write_vtu_in_parallel(filename, comm);
                      ierr = MPI_Comm_free(&comm);
                      AssertThrowMPI(ierr);
                    }
                }
              // Write in a different format than hdf5 or vtu. This case is supported, but is not
              // optimized for parallel output in that every process will write one file directly
              // into the output directory. This may or may not affect performance depending on
              // the model setup and the network file system type.
              else
                {
                  const unsigned int myid = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

                  const std::string filename = this->get_output_directory()
                                               + particles_output_base_name
                                               + "/"
                                               + particle_file_prefix
                                               + "."
                                               + Utilities::int_to_string (myid, 4)
                                               + DataOutBase::default_suffix
                                               (DataOutBase::parse_output_format(output_format));

                  std::ofstream out (filename);

                  AssertThrow(out,
                              ExcMessage("Unable to open file for writing: " + filename +"."));

                  data_out.write (out, DataOutBase::parse_output_format(output_format));
                }
            }


          // up the next time we need output
          set_last_output_time (this->get_time());

          const std::string particle_output = this->get_output_directory() + particles_output_base_name + "/" + particle_file_prefix;

          if (particle_manager == 0)
            screen_output = particle_output;

          // record the file base file name in the output file
          const std::string statistics_column_name = (particle_manager == 0 ?
                                                      "Particle file name" :
                                                      "Particle file name (" + Utilities::int_to_string(particle_manager+1) + ")");
          statistics.add_value (statistics_column_name,
                                particle_output);
        }

      return std::make_pair("Writing particle output:", screen_output);
    }



    template <int dim>
    void
    Particles<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }


    template <int dim>
    template <class Archive>
    void Particles<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &last_output_time
      & output_file_number
      & times_and_pvtu_file_names
      & output_file_names_by_timestep
      & xdmf_entries
      ;
    }


    template <int dim>
    void
    Particles<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      for (unsigned int particle_manager = 0; particle_manager < this->n_particle_managers(); ++particle_manager)
        {
          std::string particles_output_base_name = "Particles";
          if (particle_manager > 0)
            {
              particles_output_base_name += "-" + Utilities::int_to_string(particle_manager+1);
            }

          // Serialize into a stringstream. Put the following into a code
          // block of its own to ensure the destruction of the 'oa'
          // archive triggers a flush() on the stringstream so we can
          // query the completed string below.
          std::ostringstream os;
          {
            aspect::oarchive oa (os);

            this->get_particle_manager(particle_manager).save(os);
            oa << (*this);
          }

          status_strings[particles_output_base_name] = os.str();
        }
    }


    template <int dim>
    void
    Particles<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      for (unsigned int particle_manager = 0; particle_manager < this->n_particle_managers(); ++particle_manager)
        {
          std::string particles_output_base_name = "Particles";
          if (particle_manager > 0)
            {
              particles_output_base_name += "-" + Utilities::int_to_string(particle_manager+1);
            }
          // see if something was saved
          if (status_strings.find(particles_output_base_name) != status_strings.end())
            {
              std::istringstream is (status_strings.find(particles_output_base_name)->second);
              aspect::iarchive ia (is);

              // Load the particle manager
              this->get_particle_manager(particle_manager).load(is);

              ia >> (*this);
            }
        }
    }


    template <int dim>
    void
    Particles<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0.),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step.\n\n"
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");

          // now also see about the file format we're supposed to write in
          // Note: "ascii" is a legacy format used by ASPECT before particle output
          // in deal.II was implemented. It is nearly identical to the gnuplot format, thus
          // we now simply replace "ascii" by "gnuplot" should it be selected.
          prm.declare_entry ("Data output format", "vtu",
                             Patterns::MultipleSelection (DataOutBase::get_output_format_names ()+"|ascii"),
                             "A comma separated list of file formats to be used for graphical "
                             "output. The list of possible output formats that can be given "
                             "here is documented in the appendix of the manual where the current "
                             "parameter is described.");

          prm.declare_entry ("Number of grouped files", "16",
                             Patterns::Integer(0),
                             "VTU file output supports grouping files from several CPUs "
                             "into a given number of files using MPI I/O when writing on a parallel "
                             "filesystem. Select 0 for no grouping. This will disable "
                             "parallel file output and instead write one file per processor. "
                             "A value of 1 will generate one big file containing the whole "
                             "solution, while a larger value will create that many files "
                             "(at most as many as there are MPI ranks).");

          prm.declare_entry ("Write in background thread", "false",
                             Patterns::Bool(),
                             "File operations can potentially take a long time, blocking the "
                             "progress of the rest of the model run. Setting this variable to "
                             "`true' moves this process into a background thread, while the "
                             "rest of the model continues.");

          prm.declare_entry ("Temporary output location", "",
                             Patterns::Anything(),
                             "On large clusters it can be advantageous to first write the "
                             "output to a temporary file on a local file system and later "
                             "move this file to a network file system. If this variable is "
                             "set to a non-empty string it will be interpreted as a "
                             "temporary storage location.");

          prm.declare_entry ("Exclude output properties", "",
                             Patterns::Anything(),
                             "A comma separated list of particle properties that should "
                             "\\textit{not} be output. If this list contains the "
                             "entry `all', only the id of particles will be provided in "
                             "graphical output files.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      Particle::Manager<dim>::declare_parameters(prm);
    }


    template <int dim>
    void
    Particles<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particles");
        {
          output_interval = prm.get_double ("Time between data output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          AssertThrow(this->get_parameters().run_postprocessors_on_nonlinear_iterations == false,
                      ExcMessage("Postprocessing nonlinear iterations in models with "
                                 "particles is currently not supported."));

          output_formats   = Utilities::split_string_list(prm.get ("Data output format"));
          AssertThrow(Utilities::has_unique_entries(output_formats),
                      ExcMessage("The list of strings for the parameter "
                                 "'Postprocess/Particles/Data output format' contains entries more than once. "
                                 "This is not allowed. Please check your parameter file."));

          AssertThrow ((std::find (output_formats.begin(),
                                   output_formats.end(),
                                   "none") == output_formats.end())
                       ||
                       (output_formats.size() == 1),
                       ExcMessage ("If you specify 'none' for the parameter \"Data output format\", "
                                   "then this needs to be the only value given."));

          if (std::find (output_formats.begin(),
                         output_formats.end(),
                         "none") == output_formats.end())
            {
              // Note that we iterate until the value of parameters.n_particle_managers and not
              // this->particle_managers, because at this point in the program execution the
              // particle managers have not been created yet. We want to prepare as many directories
              // as there will be particle managers, once they are created.
              for (unsigned int particle_manager = 0; particle_manager < this->get_parameters().n_particle_managers; ++particle_manager)
                {
                  std::string particles_directory_base_name = "particles";
                  if (particle_manager > 0)
                    particles_directory_base_name += "-" + Utilities::int_to_string(particle_manager+1);

                  aspect::Utilities::create_directory (this->get_output_directory() + particles_directory_base_name + "/",
                                                       this->get_mpi_communicator(),
                                                       true);
                }
            }

          // Note: "ascii" is a legacy format used by ASPECT before particle output
          // in deal.II was implemented. It is nearly identical to the gnuplot format, thus
          // we simply replace "ascii" by "gnuplot" should it be selected.
          std::vector<std::string>::iterator output_format =  std::find (output_formats.begin(),
                                                                         output_formats.end(),
                                                                         "ascii");
          if (output_format != output_formats.end())
            *output_format = "gnuplot";

          group_files     = prm.get_integer("Number of grouped files");
          write_in_background_thread = prm.get_bool("Write in background thread");
          temporary_output_location = prm.get("Temporary output location");

          if (temporary_output_location != "")
            {
              // Check if a command-processor is available by calling system() with a
              // null pointer. System is guaranteed to return non-zero if it finds
              // a terminal and zero if there is none (like on the compute nodes of
              // some cluster architectures, e.g. IBM BlueGene/Q)
              AssertThrow(system((char *)nullptr) != 0,
                          ExcMessage("Usage of a temporary storage location is only supported if "
                                     "there is a terminal available to move the files to their final location "
                                     "after writing. The system() command did not succeed in finding such a terminal."));
            }

          exclude_output_properties = Utilities::split_string_list(prm.get("Exclude output properties"));

          // Never output the integrator properties that are for internal use only
          exclude_output_properties.emplace_back("internal: integrator properties");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
#define INSTANTIATE(dim) \
  template class ParticleOutput<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }

    ASPECT_REGISTER_POSTPROCESSOR(Particles,
                                  "particles",
                                  "A Postprocessor that creates particles that follow the "
                                  "velocity field of the simulation. The particles can be generated "
                                  "and propagated in various ways and they can carry a number of "
                                  "constant or time-varying properties. The postprocessor can write "
                                  "output positions and properties of all particles at chosen intervals, "
                                  "although this is not mandatory. It also allows other parts of the "
                                  "code to query the particles for information.")
  }
}
