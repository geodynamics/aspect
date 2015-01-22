/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/output.h>
#include <aspect/particle/particle.h>
#include <deal.II/numerics/data_out.h>

#ifdef DEAL_II_HAVE_HDF5
#  include <hdf5.h>
#endif

namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      // Blank class to avoid output of data, for example, if particles are used to represent
      // other physical phenomenon that influences the simulation and we don't care about their positions
      template <int dim, class T>
      class NullOutput : public Interface<dim, T>
      {
        public:
          /**
           * Constructor.
           *
           * @param[in] The directory into which output files shall be placed.
           * @param[in] The MPI communicator that describes this simulation.
           */
          NullOutput(const std::string &output_directory,
                     const MPI_Comm     communicator)
            :
            Interface<dim,T> (output_directory, communicator)
          {}

          /**
           * Write data about the particles specified in the first argument
           * to a file. If possible, encode the current simulation time
           * into this file using the data provided in the second argument.
           *
           * @param[in] particles The set of particles to generate a graphical
           *   representation for
           * @param[in] current_time Current time of the simulation, given as either
           *   years or seconds, as selected in the input file. In other words,
           *   output writers do not need to know the units in which time is
           *   described.
           * @return The name of the file that was written, or any other
           *   information that describes what output was produced if for example
           *   multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<LevelInd, T> &/*particles*/,
                               const double &/*current_time*/)
          {
            return "";
          };
      };

      template <int dim, class T>
      class ASCIIOutput : public Interface<dim, T>
      {
        public:
          /**
           * Constructor.
           *
           * @param[in] The directory into which output files shall be placed.
           * @param[in] The MPI communicator that describes this simulation.
           */
          ASCIIOutput(const std::string &output_directory,
                      const MPI_Comm     communicator)
            :
            Interface<dim,T> (output_directory, communicator)
          {}

          /**
           * Write data about the particles specified in the first argument
           * to a file. If possible, encode the current simulation time
           * into this file using the data provided in the second argument.
           *
           * @param[in] particles The set of particles to generate a graphical
           *   representation for
           * @param[in] current_time Current time of the simulation, given as either
           *   years or seconds, as selected in the input file. In other words,
           *   output writers do not need to know the units in which time is
           *   described.
           * @return The name of the file that was written, or any other
           *   information that describes what output was produced if for example
           *   multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<LevelInd, T> &particles,
                               const double &current_time)
          {
            typename std::multimap<LevelInd, T>::const_iterator  it;
            unsigned int                            i;
            std::string                             output_file_prefix, output_path_prefix, full_filename;
            std::vector<MPIDataInfo>                data_info;
            std::vector<MPIDataInfo>::iterator      dit;

            output_file_prefix = "particle-" + Utilities::int_to_string (this->file_index, 5);
            output_path_prefix = this->output_dir + output_file_prefix;
            full_filename = output_path_prefix + "." + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->communicator), 4) + ".txt";
            std::ofstream output (full_filename.c_str());
            if (!output)
              std::cout << "ERROR: proc " << Utilities::MPI::this_mpi_process(this->communicator) << " could not create " << full_filename << std::endl;

            // Get the data types
            T::add_mpi_types(data_info);

            // Print the header line
            output << "# ";
            for (dit=data_info.begin(); dit!=data_info.end(); ++dit)
              {
                // If it's a 1D element, print just the name, otherwise use []
                if (dit->n_elements == 1)
                  {
                    output << dit->name << " ";
                  }
                else
                  {
                    for (i=0; i<dit->n_elements; ++i) output << dit->name << "[" << i << "] ";
                  }
              }
            output << "\n";

            // And print the data for each particle
            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                std::vector<double>  particle_data;
                unsigned int p = 0;
                it->second.write_data(particle_data);
                for (dit=data_info.begin(); dit!=data_info.end(); ++dit)
                  {
                    for (i=0; i<dit->n_elements; ++i)
                      {
                        output << particle_data[p++] << " ";
                      }
                  }
                output << "\n";
              }

            output.close();

            this->file_index++;

            return output_path_prefix;
          };
      };



      template <int dim, class T>
      class VTUOutput : public Interface<dim, T>
      {
        private:

          /**
           * A list of pairs (time, pvtu_filename) that have so far been written
           * and that we will pass to DataOutInterface::write_pvd_record
           * to create a master file that can make the association
           * between simulation time and corresponding file name (this
           * is done because there is no way to store the simulation
           * time inside the .pvtu or .vtu files).
           */
//TODO: This needs to be serialized
          std::vector<std::pair<double,std::string> > times_and_pvtu_file_names;

          /**
           * Like the previous variable, but for the .visit file.
           */
          std::vector<std::string>                    vtu_file_names;

        public:
          /**
           * Constructor.
           *
           * @param[in] The directory into which output files shall be placed.
           * @param[in] The MPI communicator that describes this simulation.
           */
          VTUOutput(const std::string &output_directory,
                    const MPI_Comm     communicator)
            :
            Interface<dim,T> (output_directory, communicator)
          {}

          /**
           * Write data about the particles specified in the first argument
           * to a file. If possible, encode the current simulation time
           * into this file using the data provided in the second argument.
           *
           * @param[in] particles The set of particles to generate a graphical
           *   representation for
           * @param[in] current_time Current time of the simulation, given as either
           *   years or seconds, as selected in the input file. In other words,
           *   output writers do not need to know the units in which time is
           *   described.
           * @return The name of the file that was written, or any other
           *   information that describes what output was produced if for example
           *   multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<LevelInd, T> &particles,
                               const double &current_time)
          {
            std::vector<MPIDataInfo>                data_info;
            unsigned int                            data_offset;

            const std::string output_file_prefix = "particles-" + Utilities::int_to_string (this->file_index, 5);
            const std::string output_path_prefix = this->output_dir + output_file_prefix;

            const std::string filename = (output_file_prefix +
                                          "." +
                                          Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->communicator), 4) +
                                          ".vtu");
            const std::string full_filename = (this->output_dir + filename);

            std::ofstream output (full_filename.c_str());
            AssertThrow (output, ExcIO());

            const unsigned int n_particles = particles.size();

            // Write VTU file XML
            output << "<?xml version=\"1.0\"?>\n";
            output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            output << "  <UnstructuredGrid>\n";
            output << "    <Piece NumberOfPoints=\"" << n_particles << "\" NumberOfCells=\"" << n_particles << "\">\n";

            // Go through the particles on this domain and print the position of each one
            output << "      <Points>\n";
            output << "        <DataArray name=\"Position\" type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
            for (typename std::multimap<LevelInd, T>::const_iterator
                 it=particles.begin(); it!=particles.end(); ++it)
              {
                output << "          " << it->second.get_location();

                // pad with zeros since VTU format wants x/y/z coordinates
                for (unsigned int d=dim; d<3; ++d)
                  output << " 0.0";

                output << "\n";
              }
            output << "        </DataArray>\n";
            output << "      </Points>\n";

            // Write cell related data (empty)
            output << "      <Cells>\n";
            output << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
            for (unsigned int i=1; i<=n_particles; ++i)
              output << "          " << i-1 << "\n";
            output << "        </DataArray>\n";
            output << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
            for (unsigned int i=1; i<=n_particles; ++i)
              output << "          " << i << "\n";
            output << "        </DataArray>\n";
            output << "        <DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">\n";
            for (unsigned int i=1; i<=n_particles; ++i)
              output << "          1\n";
            output << "        </DataArray>\n";
            output << "      </Cells>\n";

            // Get the data types
            T::add_mpi_types(data_info);

            // Write data for each particle (id, velocity, etc)
            output << "      <PointData Scalars=\"scalars\">\n";

            // Print the data associated with the particles, skipping the first entry (position)

            std::vector<MPIDataInfo>::iterator dit = data_info.begin();
            data_offset = dit->n_elements;
            dit++;
            for (; dit!=data_info.end(); ++dit)
              {
                output << "        <DataArray type=\"Float64\" Name=\"" << dit->name << "\" NumberOfComponents=\"" << (dit->n_elements == 2 ? 3 : dit->n_elements) << "\" Format=\"ascii\">\n";
                for (typename std::multimap<LevelInd, T>::const_iterator
                     it=particles.begin(); it!=particles.end(); ++it)
                  {
                    std::vector<double> particle_data;
                    it->second.write_data(particle_data);
                    output << "          ";
                    for (unsigned int d=0; d<dit->n_elements; ++d)
                      {
                        output << particle_data[data_offset+d] << " ";
                      }
                    if (dit->n_elements == 2)
                      output << "0 ";
                    output << "\n";
                  }
                data_offset += dit->n_elements;
                output << "        </DataArray>\n";
              }
            output << "      </PointData>\n";

            output << "    </Piece>\n";
            output << "  </UnstructuredGrid>\n";
            output << "</VTKFile>\n";

            output.close();


            // Write the parallel pvtu and pvd files on the root process
            if (Utilities::MPI::this_mpi_process(this->communicator) == 0)
              {
                const std::string pvtu_filename = (output_file_prefix + ".pvtu");
                const std::string full_pvtu_filename = (this->output_dir + pvtu_filename);

                std::ofstream pvtu_output (full_pvtu_filename.c_str());
                AssertThrow (pvtu_output, ExcIO());

                pvtu_output << "<?xml version=\"1.0\"?>\n";
                pvtu_output << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                pvtu_output << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
                pvtu_output << "    <PPoints>\n";
                pvtu_output << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
                pvtu_output << "    </PPoints>\n";
                pvtu_output << "    <PPointData Scalars=\"scalars\">\n";

                for (dit=data_info.begin()+1; dit!=data_info.end(); ++dit)
                  {
                    pvtu_output << "      <PDataArray type=\"Float64\" Name=\"" << dit->name << "\" NumberOfComponents=\"" << (dit->n_elements == 2 ? 3 : dit->n_elements) << "\" format=\"ascii\"/>\n";
                  }
                pvtu_output << "    </PPointData>\n";
                for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(this->communicator); ++i)
                  {
                    pvtu_output << "    <Piece Source=\"" << output_file_prefix << "." << Utilities::int_to_string(i, 4) << ".vtu\"/>\n";
                  }
                pvtu_output << "  </PUnstructuredGrid>\n";
                pvtu_output << "</VTKFile>\n";
                pvtu_output.close();

                times_and_pvtu_file_names.push_back(std::make_pair(current_time,
                                                                   pvtu_filename));
                vtu_file_names.push_back (full_filename);

                // write .pvd and .visit records
                const std::string pvd_master_filename = (this->output_dir + "particles.pvd");
                std::ofstream pvd_master (pvd_master_filename.c_str());
                DataOut<dim>().write_pvd_record (pvd_master, times_and_pvtu_file_names);

//TODO: write a global .visit record. this needs a variant of the write_visit_record
// function
                /*
                              const std::string visit_master_filename = (this->output_dir + "particles.visit");
                              std::ofstream visit_master (visit_master_filename.c_str());
                              DataOut<dim>().write_visit_record (visit_master, vtu_file_names);
                */
              }
            this->file_index++;

            return output_path_prefix;
          }
      };



      template <int dim, class T>
      class HDF5Output : public Interface<dim, T>
      {
        private:
          // A set of XDMF data objects to create the XDMF file for particles
          std::vector<XDMFEntry>          xdmf_entries;

        public:
          /**
           * Constructor.
           *
           * @param[in] The directory into which output files shall be placed.
           * @param[in] The MPI communicator that describes this simulation.
           */
          HDF5Output(const std::string &output_directory,
                     const MPI_Comm     communicator)
            :
            Interface<dim,T> (output_directory, communicator)
          {}

          /**
           * Write data about the particles specified in the first argument
           * to a file. If possible, encode the current simulation time
           * into this file using the data provided in the second argument.
           *
           * @param[in] particles The set of particles to generate a graphical
           *   representation for
           * @param[in] current_time Current time of the simulation, given as either
           *   years or seconds, as selected in the input file. In other words,
           *   output writers do not need to know the units in which time is
           *   described.
           * @return The name of the file that was written, or any other
           *   information that describes what output was produced if for example
           *   multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<LevelInd, T> &particles,
                               const double &current_time)
          {
            std::string             output_file_prefix, output_path_prefix, full_filename;

            output_file_prefix = "particle-" + Utilities::int_to_string (this->file_index, 5);
            output_path_prefix = this->output_dir + output_file_prefix;
#ifdef DEAL_II_HAVE_HDF5
            // TODO: error checking for H5 calls
            hid_t   h5_file_id;
            hid_t   plist_id, xfer_plist_id, position_dataset_id, velocity_dataset_id, pid_dataset_id;
            hid_t   file_dataspace_id, pos_file_dataspace_id, vel_file_dataspace_id, pid_file_dataspace_id;
            hid_t   dim_dataspace_id, one_dim_ds_id;
            hid_t   dim_mem_ds_id, one_dim_mem_ds_id;
            hid_t   pattr_id;
            hsize_t dims[2], offset[2], count[2];
            unsigned int  mpi_offset, mpi_count, local_particle_count, global_particle_count, d, i;

            double  *pos_data, *vel_data, *id_data;

            std::string h5_filename = output_path_prefix+".h5";

            // Create parallel file access
            plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef H5_HAVE_PARALLEL
            H5Pset_fapl_mpio(plist_id, this->communicator, MPI_INFO_NULL);
#endif

            // Create the file
            h5_file_id = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

            // Create the file dataspace descriptions
            local_particle_count = particles.size();
            // TODO: error checking on MPI call
            MPI_Comm com = this->communicator;
            MPI_Allreduce(&local_particle_count, &global_particle_count, 1, MPI_UNSIGNED, MPI_SUM, com);
            dims[0] = global_particle_count;
            dims[1] = 3;
            one_dim_ds_id = H5Screate_simple(1, dims, NULL);
            dim_dataspace_id = H5Screate_simple(2, dims, NULL);

            // Create the datasets
#if H5Dcreate_vers == 1
            position_dataset_id = H5Dcreate(h5_file_id, "nodes", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT);
            velocity_dataset_id = H5Dcreate(h5_file_id, "velocity", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT);
            pid_dataset_id = H5Dcreate(h5_file_id, "id", H5T_NATIVE_DOUBLE, one_dim_ds_id, H5P_DEFAULT);
#else
            position_dataset_id = H5Dcreate(h5_file_id, "nodes", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            velocity_dataset_id = H5Dcreate(h5_file_id, "velocity", H5T_NATIVE_DOUBLE, dim_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            pid_dataset_id = H5Dcreate(h5_file_id, "id", H5T_NATIVE_DOUBLE, one_dim_ds_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

            // Close the file dataspaces
            H5Sclose(dim_dataspace_id);
            H5Sclose(one_dim_ds_id);

            // Just in case we forget what they are
            count[0] = 1;
            file_dataspace_id = H5Screate_simple (1, count, NULL);
#if H5Acreate_vers == 1
            pattr_id = H5Acreate(h5_file_id, "Ermahgerd! Pertecrs!", H5T_NATIVE_INT, file_dataspace_id, H5P_DEFAULT);
#else
            pattr_id = H5Acreate(h5_file_id, "Ermahgerd! Pertecrs!", H5T_NATIVE_INT, file_dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
#endif
            H5Aclose(pattr_id);
            H5Sclose(file_dataspace_id);

            // Count the number of particles on this process and get the offset among all processes
            mpi_count = particles.size();
            MPI_Scan(&mpi_count, &mpi_offset, 1, MPI_UNSIGNED, MPI_SUM, this->communicator);
            count[0] = mpi_count;
            count[1] = 3;
            offset[0] = mpi_offset-mpi_count;
            offset[1] = 0;

            // Select the appropriate dataspace for this process
            one_dim_mem_ds_id = H5Screate_simple(1, count, NULL);
            dim_mem_ds_id = H5Screate_simple(2, count, NULL);
            pos_file_dataspace_id = H5Dget_space(position_dataset_id);
            vel_file_dataspace_id = H5Dget_space(velocity_dataset_id);
            pid_file_dataspace_id = H5Dget_space(pid_dataset_id);

            // And select the hyperslabs from each dataspace
            H5Sselect_hyperslab(pos_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(vel_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(pid_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            // Create property list for collective dataset write
            xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
            H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
#endif

            // Read the local particle data
            pos_data = new double[3*particles.size()];
            vel_data = new double[3*particles.size()];
            id_data = new double[particles.size()];

            typename std::multimap<LevelInd, T>::const_iterator it;
            for (i=0,it=particles.begin(); it!=particles.end(); ++i,++it)
              {
                for (d=0; d<dim; ++d)
                  {
                    pos_data[i*3+d] = it->second.get_location()(d);
                    vel_data[i*3+d] = it->second.get_velocity()(d);
                  }
                if (dim < 3)
                  {
                    pos_data[i*3+2] = 0;
                    vel_data[i*3+2] = 0;
                  }
                id_data[i] = it->second.get_id();
              }

            // Write particle data to the HDF5 file
            H5Dwrite(position_dataset_id, H5T_NATIVE_DOUBLE, dim_mem_ds_id, pos_file_dataspace_id, xfer_plist_id, pos_data);
            H5Dwrite(velocity_dataset_id, H5T_NATIVE_DOUBLE, dim_mem_ds_id, vel_file_dataspace_id, xfer_plist_id, vel_data);
            H5Dwrite(pid_dataset_id, H5T_NATIVE_DOUBLE, one_dim_mem_ds_id, pid_file_dataspace_id, xfer_plist_id, id_data);

            // Clear allocated resources
            delete[] pos_data;
            delete[] vel_data;
            delete[] id_data;

            H5Pclose(xfer_plist_id);
            H5Sclose(one_dim_mem_ds_id);
            H5Sclose(dim_mem_ds_id);
            H5Sclose(pos_file_dataspace_id);
            H5Sclose(vel_file_dataspace_id);
            H5Sclose(pid_file_dataspace_id);
            H5Dclose(position_dataset_id);
            H5Dclose(velocity_dataset_id);
            H5Dclose(pid_dataset_id);
            H5Pclose(plist_id);
            H5Fclose(h5_file_id);

            // Record and output XDMF info on root process
            if (Utilities::MPI::this_mpi_process(this->communicator) == 0)
              {
                std::string local_h5_filename = output_file_prefix+".h5";
                XDMFEntry   entry(local_h5_filename, current_time, global_particle_count, 0, 3);
                DataOut<dim> data_out;
                const std::string xdmf_filename = (this->output_dir + "particle.xdmf");

                entry.add_attribute("velocity", 3);
                entry.add_attribute("id", 1);
                xdmf_entries.push_back(entry);

                data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), this->communicator);
              }
#endif

            this->file_index++;

            return output_path_prefix;
          }
      };



      template <int dim, class T>
      Interface<dim,T> *
      create_output_object (const std::string &data_output_format,
                            const std::string &output_directory,
                            const MPI_Comm     communicator)
      {
        if (data_output_format == "ascii")
          return new ASCIIOutput<dim,T>(output_directory, communicator);
        else if (data_output_format == "vtu")
          return new VTUOutput<dim,T>(output_directory, communicator);
        else if (data_output_format == "hdf5")
          return new HDF5Output<dim,T>(output_directory, communicator);
        else if (data_output_format == "none")
          return new NullOutput<dim,T>(output_directory, communicator);
        else
          Assert (false, ExcNotImplemented());

        return 0;
      }


      std::string
      output_object_names ()
      {
        return ("none|"
                "ascii|"
                "vtu|"
                "hdf5");
      }


      // explicit instantiations
      template
      Interface<2,Particle::BaseParticle<2> > *
      create_output_object (const std::string &data_output_format,
                            const std::string &output_directory,
                            const MPI_Comm     communicator);
      template
      Interface<3,Particle::BaseParticle<3> > *
      create_output_object (const std::string &data_output_format,
                            const std::string &output_directory,
                            const MPI_Comm     communicator);
    }
  }
}
