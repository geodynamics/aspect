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

#include <aspect/particle/output/hdf5.h>

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
      template <int dim>
      HDF5Output<dim>::HDF5Output()
        :
        file_index(0)
      {}

      template <int dim>
      std::string
      HDF5Output<dim>::output_particle_data(const std::multimap<LevelInd, Particle<dim> > &particles,
                                            const std::vector<std::string>  &/*data_names*/,
                                            const std::vector<unsigned int> &/*data_components*/,
                                            const double &current_time)
      {
        // avoid warnings about unused variables
        (void)current_time;

        std::string             output_file_prefix, output_path_prefix, full_filename;

        output_file_prefix = "particle-" + Utilities::int_to_string (file_index, 5);
        output_path_prefix = this->get_output_directory() + output_file_prefix;
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
        H5Pset_fapl_mpio(plist_id, this->get_mpi_communicator(), MPI_INFO_NULL);
#endif

        // Create the file
        h5_file_id = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

        // Create the file dataspace descriptions
        local_particle_count = particles.size();
        // TODO: error checking on MPI call
        MPI_Comm com = this->get_mpi_communicator();
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
        MPI_Scan(&mpi_count, &mpi_offset, 1, MPI_UNSIGNED, MPI_SUM, this->get_mpi_communicator());
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

        typename std::multimap<LevelInd, Particle<dim> >::const_iterator it;
        for (i=0,it=particles.begin(); it!=particles.end(); ++i,++it)
          {
            for (d=0; d<dim; ++d)
              {
                pos_data[i*3+d] = it->second.get_location()(d);
                vel_data[i*3+d] = 0; //it->second.get_velocity()(d);
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
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            std::string local_h5_filename = output_file_prefix+".h5";
            XDMFEntry   entry(local_h5_filename, current_time, global_particle_count, 0, 3);
            DataOut<dim> data_out;
            const std::string xdmf_filename = (this->get_output_directory() + "particle.xdmf");

            entry.add_attribute("velocity", 3);
            entry.add_attribute("id", 1);
            xdmf_entries.push_back(entry);

            data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), this->get_mpi_communicator());
          }

#else
        // avoid warnings about unused variables
        (void)particles;
        (void)current_time;
#endif

        file_index++;

        return output_path_prefix;
      }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      ASPECT_REGISTER_PARTICLE_OUTPUT(HDF5Output,
                                      "hdf5",
                                      "This particle output plugin writes particle "
                                      "positions and properties into hdf5 files.")
    }
  }
}

