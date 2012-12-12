/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */

#ifndef __aspect__particle_output_h
#define __aspect__particle_output_h

#include <deal.II/numerics/data_out.h>
#include <aspect/particle/particle.h>

#ifdef DEAL_II_HAVE_HDF5
#include <hdf5.h>
#endif

namespace aspect
{
  namespace Particle
  {
    // Abstract class representing particle output
    template <int dim, class T>
    class Output
    {
      protected:
        // MPI communicator to be used for output synchronization
        MPI_Comm        _communicator;

        // Internal index of file output number, must be incremented by subclasses when they create a new file
        unsigned int    _file_index;

        // Path to directory in which to put particle output files
        std::string     _output_dir;

      public:
        Output(void)
        {
          _file_index = 0;
        };
        virtual ~Output(void) {};

        unsigned int self_rank(void)
        {
          return Utilities::MPI::this_mpi_process(_communicator);
        };

        unsigned int world_size(void)
        {
          return Utilities::MPI::n_mpi_processes(_communicator);
        };

        void set_mpi_comm(MPI_Comm new_comm_world)
        {
          _communicator = new_comm_world;
        };

        void set_output_directory(std::string new_out_dir)
        {
          _output_dir = new_out_dir;
        };
        virtual std::string output_particle_data(const std::multimap<LevelInd, T> &particles, const double &current_time) = 0;
    };

    template <int dim, class T>
    class NullOutput : public Output<dim, T>
    {
      public:
        virtual std::string output_particle_data(const std::multimap<LevelInd, T> &particles, const double &current_time)
        {
          return "";
        };
    };

    template <int dim, class T>
    class ASCIIOutput : public Output<dim, T>
    {
      public:
        virtual std::string output_particle_data(const std::multimap<LevelInd, T> &particles, const double &current_time)
        {
          typename std::multimap<LevelInd, T>::const_iterator  it;
          unsigned int                            i;
          std::string                             output_file_prefix, output_path_prefix, full_filename;
          std::vector<MPIDataInfo>                data_info;
          std::vector<MPIDataInfo>::iterator      dit;
          char                                    *particle_data, *p;

          output_file_prefix = "particle-" + Utilities::int_to_string (Output<dim, T>::_file_index, 5);
          output_path_prefix = Output<dim, T>::_output_dir + output_file_prefix;
          full_filename = output_path_prefix + "." + Utilities::int_to_string(Output<dim, T>::self_rank(), 4) + ".txt";
          std::ofstream output (full_filename.c_str());
          if (!output) std::cout << "ERROR: proc " << Output<dim, T>::self_rank() << " could not create " << full_filename << std::endl;

          // Get the data types
          T::add_mpi_types(data_info);

          // Print the header line
          output << "# ";
          for (dit=data_info.begin(); dit!=data_info.end(); ++dit)
            {
              // If it's a 1D element, print just the name, otherwise use []
              if (dit->_num_elems == 1)
                {
                  output << dit->_name << " ";
                }
              else
                {
                  for (i=0; i<dit->_num_elems; ++i) output << dit->_name << "[" << i << "] ";
                }
            }
          output << "\n";

          // And print the data for each particle
          particle_data = new char[T::data_len(HDF5_DATA)];
          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              it->second.write_data(HDF5_DATA, particle_data);
              p = particle_data;
              for (dit=data_info.begin(); dit!=data_info.end(); ++dit)
                {
                  // Currently assumes all data is double, may need to change this
                  for (i=0; i<dit->_num_elems; ++i)
                    {
                      output << ((double *)p)[0] << " ";
                      p += sizeof(double);
                    }
                }
              output << "\n";
            }
          delete[] particle_data;

          output.close();

          Output<dim, T>::_file_index++;

          return output_path_prefix;
        };
    };

    template <int dim, class T>
    class VTUOutput : public Output<dim, T>
    {
        /**
         * A list of pairs (time, pvtu_filename) that have so far been written
         * and that we will pass to DataOutInterface::write_pvd_record
         * to create a master file that can make the association
         * between simulation time and corresponding file name (this
         * is done because there is no way to store the simulation
         * time inside the .pvtu or .vtu files).
         */
        std::vector<std::pair<double,std::string> > times_and_pvtu_names;

      public:
        virtual std::string output_particle_data(const std::multimap<LevelInd, T> &particles, const double &current_time)
        {
          typename std::multimap<LevelInd, T>::const_iterator  it;
          unsigned int                            d, i, num_particles;
          DataOut<dim>                            data_out;
          std::vector<MPIDataInfo>                data_info;
          std::vector<MPIDataInfo>::iterator      dit;
          std::string                             output_file_prefix, output_path_prefix;
          unsigned int                            data_offset;
          char                                    *particle_data, *p;

          output_file_prefix = "particle-" + Utilities::int_to_string (Output<dim, T>::_file_index, 5);
          output_path_prefix = Output<dim, T>::_output_dir + output_file_prefix;

          const std::string filename = (output_file_prefix +
                                        "." +
                                        Utilities::int_to_string(Output<dim, T>::self_rank(), 4) +
                                        ".vtu");
          const std::string full_filename = (Output<dim, T>::_output_dir + filename);

          std::ofstream output (full_filename.c_str());
          if (!output) std::cout << "ERROR: proc " << Output<dim, T>::self_rank() << " could not create " << filename << std::endl;

          num_particles = particles.size();

          // Write VTU file XML
          output << "<?xml version=\"1.0\"?>\n";
          output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
          output << "  <UnstructuredGrid>\n";
          output << "    <Piece NumberOfPoints=\"" << num_particles << "\" NumberOfCells=\"" << num_particles << "\">\n";

          // Go through the particles on this domain and print the position of each one
          // TODO: can we change this to dim-dimensions rather than 3?
          output << "      <Points>\n";
          output << "        <DataArray name=\"Position\" type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              output << "          ";
              for (d=0; d<3; ++d)
                {
                  if (d < dim) output << it->second.location()[d] << " ";
                  else output << "0.0 ";
                }
              output << "\n";
            }
          output << "        </DataArray>\n";
          output << "      </Points>\n";

          // Write cell related data (empty)
          output << "      <Cells>\n";
          output << "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";
          for (i=1; i<=num_particles; ++i) output << "          " << i-1 << "\n";
          output << "        </DataArray>\n";
          output << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
          for (i=1; i<=num_particles; ++i) output << "          " << i << "\n";
          output << "        </DataArray>\n";
          output << "        <DataArray type=\"UInt8\" Name=\"types\" Format=\"ascii\">\n";
          for (i=1; i<=num_particles; ++i) output << "          1\n";
          output << "        </DataArray>\n";
          output << "      </Cells>\n";

          // Get the data types
          T::add_mpi_types(data_info);

          // Write data for each particle (id, velocity, etc)
          output << "      <PointData Scalars=\"scalars\">\n";

          // Print the data associated with the particles, skipping the first entry (assumed to be position)
          particle_data = new char[T::data_len(HDF5_DATA)];
          dit = data_info.begin();
          data_offset = dit->_num_elems*dit->_elem_size_bytes;
          dit++;
          for (; dit!=data_info.end(); ++dit)
            {
              output << "        <DataArray type=\"Float64\" Name=\"" << dit->_name << "\" NumberOfComponents=\"" << (dit->_num_elems == 2 ? 3 : dit->_num_elems) << "\" Format=\"ascii\">\n";
              for (it=particles.begin(); it!=particles.end(); ++it)
                {
                  it->second.write_data(HDF5_DATA, particle_data);
                  p = particle_data+data_offset;
                  output << "          ";
                  for (d=0; d<dit->_num_elems; ++d)
                    {
                      output << ((double *)p)[0] << " ";
                      p += sizeof(double);
                    }
                  if (dit->_num_elems == 2) output << "0 ";
                  output << "\n";
                }
              data_offset += dit->_num_elems*dit->_elem_size_bytes;
              output << "        </DataArray>\n";
            }
          delete[] particle_data;
          output << "      </PointData>\n";

          output << "    </Piece>\n";
          output << "  </UnstructuredGrid>\n";
          output << "</VTKFile>\n";

          output.close();

          // Write the parallel pvtu and pvd files on the root process
          if (Output<dim, T>::self_rank() == 0)
            {
              const std::string pvtu_filename = (output_path_prefix + ".pvtu");

              std::ofstream pvtu_output (pvtu_filename.c_str());
              if (!pvtu_output) std::cout << "ERROR: Could not create " << filename << std::endl;

              pvtu_output << "<?xml version=\"1.0\"?>\n";
              pvtu_output << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
              pvtu_output << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
              pvtu_output << "    <PPoints>\n";
              pvtu_output << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
              pvtu_output << "    </PPoints>\n";
              pvtu_output << "    <PPointData Scalars=\"scalars\">\n";

              for (dit=data_info.begin()+1; dit!=data_info.end(); ++dit)
                {
                  pvtu_output << "      <PDataArray type=\"Float64\" Name=\"" << dit->_name << "\" NumberOfComponents=\"" << (dit->_num_elems == 2 ? 3 : dit->_num_elems) << "\" format=\"ascii\"/>\n";
                }
              pvtu_output << "    </PPointData>\n";
              for (i=0; i<Output<dim, T>::world_size(); ++i)
                {
                  pvtu_output << "    <Piece Source=\"" << output_file_prefix << "." << Utilities::int_to_string(i, 4) << ".vtu\"/>\n";
                }
              pvtu_output << "  </PUnstructuredGrid>\n";
              pvtu_output << "</VTKFile>\n";
              pvtu_output.close();

              times_and_pvtu_names.push_back(std::pair<double,std::string>(current_time, output_file_prefix+".pvtu"));
              const std::string pvd_master_filename = (Output<dim, T>::_output_dir + "particle.pvd");
              std::ofstream pvd_master (pvd_master_filename.c_str());
              data_out.write_pvd_record (pvd_master, times_and_pvtu_names);
              pvd_master.close();
            }
          Output<dim, T>::_file_index++;

          return output_path_prefix;
        };
    };

    template <int dim, class T>
    class HDF5Output : public Output<dim, T>
    {
        // A set of XDMF data objects to create the XDMF file for particles
        std::vector<XDMFEntry>          xdmf_entries;

      public:
        virtual std::string output_particle_data(const std::multimap<LevelInd, T> &particles, const double &current_time)
        {
          typename std::multimap<LevelInd, T>::const_iterator  it;
          unsigned int            d, i;
          std::string             output_file_prefix, output_path_prefix, full_filename;

          output_file_prefix = "particle-" + Utilities::int_to_string (Output<dim, T>::_file_index, 5);
          output_path_prefix = Output<dim, T>::_output_dir + output_file_prefix;
#ifdef DEAL_II_HAVE_HDF5
          // TODO: error checking for H5 calls
          hid_t   h5_file_id;
          hid_t   plist_id, xfer_plist_id, position_dataset_id, velocity_dataset_id, pid_dataset_id;
          hid_t   file_dataspace_id, pos_file_dataspace_id, vel_file_dataspace_id, pid_file_dataspace_id;
          hid_t   dim_dataspace_id, one_dim_ds_id;
          hid_t   dim_mem_ds_id, one_dim_mem_ds_id;
          hid_t   pattr_id;
          herr_t  status;
          hsize_t dims[2], offset[2], count[2];
          unsigned int  mpi_offset, mpi_count, local_particle_count, global_particle_count;

          double  *pos_data, *vel_data, *id_data;

          std::string h5_filename = output_path_prefix+".h5";

          // Create parallel file access
          plist_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef H5_HAVE_PARALLEL
          H5Pset_fapl_mpio(plist_id, Output<dim, T>::_communicator, MPI_INFO_NULL);
#endif

          // Create the file
          h5_file_id = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

          // Create the file dataspace descriptions
          local_particle_count = particles.size();
          // TODO: error checking on MPI call
          MPI_Allreduce(&local_particle_count, &global_particle_count, 1, MPI_UNSIGNED, MPI_SUM, Output<dim, T>::_communicator);
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
          MPI_Scan(&mpi_count, &mpi_offset, 1, MPI_UNSIGNED, MPI_SUM, Output<dim, T>::_communicator);
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
          status = H5Sselect_hyperslab(pos_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          status = H5Sselect_hyperslab(vel_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          status = H5Sselect_hyperslab(pid_file_dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

          // Create property list for collective dataset write
          xfer_plist_id = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
          H5Pset_dxpl_mpio(xfer_plist_id, H5FD_MPIO_COLLECTIVE);
#endif

          // Read the local particle data
          pos_data = new double[3*particles.size()];
          vel_data = new double[3*particles.size()];
          id_data = new double[particles.size()];

          for (i=0,it=particles.begin(); it!=particles.end(); ++i,++it)
            {
              for (d=0; d<dim; ++d)
                {
                  pos_data[i*3+d] = it->second.location()(d);
                  vel_data[i*3+d] = it->second.velocity()(d);
                }
              if (dim < 3)
                {
                  pos_data[i*3+2] = 0;
                  vel_data[i*3+2] = 0;
                }
              id_data[i] = it->second.id_num();
            }

          // Write particle data to the HDF5 file
          H5Dwrite(position_dataset_id, H5T_NATIVE_DOUBLE, dim_mem_ds_id, pos_file_dataspace_id, xfer_plist_id, pos_data);
          H5Dwrite(velocity_dataset_id, H5T_NATIVE_DOUBLE, dim_mem_ds_id, vel_file_dataspace_id, xfer_plist_id, vel_data);
          H5Dwrite(pid_dataset_id, H5T_NATIVE_DOUBLE, one_dim_mem_ds_id, pid_file_dataspace_id, xfer_plist_id, id_data);

          // Clear allocated resources
          delete[] pos_data;
          delete[] vel_data;
          delete[] id_data;
          status = H5Pclose(xfer_plist_id);
          status = H5Sclose(one_dim_mem_ds_id);
          status = H5Sclose(dim_mem_ds_id);
          status = H5Sclose(pos_file_dataspace_id);
          status = H5Sclose(vel_file_dataspace_id);
          status = H5Sclose(pid_file_dataspace_id);
          status = H5Dclose(position_dataset_id);
          status = H5Dclose(velocity_dataset_id);
          status = H5Dclose(pid_dataset_id);
          status = H5Pclose(plist_id);
          status = H5Fclose(h5_file_id);

          // Record and output XDMF info on root process
          if (Output<dim, T>::self_rank() == 0)
            {
              std::string local_h5_filename = output_file_prefix+".h5";
              XDMFEntry   entry(local_h5_filename, current_time, global_particle_count, 0, 3);
              DataOut<dim> data_out;
              const std::string xdmf_filename = (Output<dim, T>::_output_dir + "particle.xdmf");

              entry.add_attribute("velocity", 3);
              entry.add_attribute("id", 1);
              xdmf_entries.push_back(entry);

              data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), Output<dim, T>::_communicator);
            }
#endif

          Output<dim, T>::_file_index++;

          return output_path_prefix;
        };
    };
  }
}

#endif
