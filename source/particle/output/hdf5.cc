/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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

#ifdef DEAL_II_WITH_HDF5
#include <hdf5.h>
#endif

namespace aspect
{
  namespace Particle
  {
    namespace Output
    {

      // Define the hdf5 type of the partice index output
#ifdef DEAL_II_WITH_HDF5

#ifdef DEAL_II_WITH_64BIT_INDICES
#  define HDF5_PARTICLE_INDEX_TYPE H5T_NATIVE_ULLONG
#else
#  define HDF5_PARTICLE_INDEX_TYPE H5T_NATIVE_UINT
#endif

#endif

      template <int dim>
      HDF5Output<dim>::HDF5Output()
        :
        file_index(0)
      {

#ifndef DEAL_II_WITH_HDF5
        AssertThrow (false,
                     ExcMessage ("deal.ii was not compiled with HDF5 support, "
                                 "so HDF5 output is not possible. Please "
                                 "recompile deal.ii with HDF5 support turned on "
                                 "or select a different tracer output format."));
#endif

      }

      template <int dim>
      std::string
      HDF5Output<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                            const std::vector<std::pair<std::string, unsigned int> > &property_component_list,
                                            const double current_time)
      {
#ifdef DEAL_II_WITH_HDF5
        // Create the filename
        const std::string output_file_prefix = "particle-" + Utilities::int_to_string (file_index, 5);
        const std::string output_path_prefix = this->get_output_directory() + output_file_prefix;
        const std::string h5_filename = output_path_prefix+".h5";

        // Create the hdf5 output size information
        types::particle_index n_local_particles = particles.size();
        const types::particle_index n_global_particles = Utilities::MPI::sum(n_local_particles,this->get_mpi_communicator());

        hsize_t global_dataset_size[2];
        global_dataset_size[0] = n_global_particles;
        global_dataset_size[1] = 3;

        hsize_t local_dataset_size[2];
        local_dataset_size[0] = n_local_particles;
        local_dataset_size[1] = 3;

        // Get the offset of the local particles among all processes
        types::particle_index local_particle_index_offset;
        MPI_Scan(&n_local_particles, &local_particle_index_offset, 1, ASPECT_TRACER_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());

        hsize_t offset[2];
        offset[0] = local_particle_index_offset - n_local_particles;
        offset[1] = 0;

        // Prepare the output data
        std::vector<double> position_data (3 * n_local_particles,0.0);
        std::vector<types::particle_index> index_data (n_local_particles);
        std::vector<std::vector<double> > property_data(property_component_list.size());

        for (unsigned int property = 0; property < property_component_list.size(); ++property)
          {
            const unsigned int data_components = (property_component_list[property].second != 2
                                                  ?
                                                  property_component_list[property].second
                                                  :
                                                  3);

            property_data[property].resize(data_components * n_local_particles,0.0);
          }

        typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator it = particles.begin();
        for (unsigned int i = 0; it != particles.end(); ++i, ++it)
          {
            for (unsigned int d = 0; d < dim; ++d)
              position_data[i*3+d] = it->second.get_location()(d);

            index_data[i] = it->second.get_id();

            const std::vector<double> properties = it->second.get_properties();
            unsigned int particle_property_index = 0;

            for (unsigned int property = 0; property < property_component_list.size(); ++property)
              {

                const unsigned int data_components = (property_component_list[property].second != 2
                                                      ?
                                                      property_component_list[property].second
                                                      :
                                                      3);

                for (unsigned int component = 0; component < property_component_list[property].second; ++component,++particle_property_index)
                  property_data[property][i * data_components + component] = properties[particle_property_index];

              }
          }

        // Create parallel file access
        const hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        // Create property list for collective dataset write
        const hid_t write_properties = H5Pcreate(H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
        H5Pset_fapl_mpio(plist_id, this->get_mpi_communicator(), MPI_INFO_NULL);
        H5Pset_dxpl_mpio(write_properties, H5FD_MPIO_COLLECTIVE);
#endif

        // Create the file
        const hid_t h5_file = H5Fcreate(h5_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        H5Pclose(plist_id);

        // Write the position data
        const hid_t position_dataspace = H5Screate_simple(2, global_dataset_size, NULL);
        const hid_t local_position_dataspace = H5Screate_simple(2, local_dataset_size, NULL);

#if H5Dcreate_vers == 1
        hid_t position_dataset = H5Dcreate(h5_file, "nodes", H5T_NATIVE_DOUBLE, position_dataspace, H5P_DEFAULT);
#else
        hid_t position_dataset = H5Dcreate(h5_file, "nodes", H5T_NATIVE_DOUBLE, position_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

        // Select the local hyperslab from the dataspace
        H5Sselect_hyperslab(position_dataspace, H5S_SELECT_SET, offset, NULL, local_dataset_size, NULL);

        // Write position data to the HDF5 file
        H5Dwrite(position_dataset, H5T_NATIVE_DOUBLE, local_position_dataspace, position_dataspace, write_properties, &position_data[0]);

        H5Sclose(local_position_dataspace);
        H5Sclose(position_dataspace);
        H5Dclose(position_dataset);

        // Write the index data
        const hid_t index_dataspace = H5Screate_simple(1, global_dataset_size, NULL);
        const hid_t local_index_dataspace = H5Screate_simple(1, local_dataset_size, NULL);
#if H5Dcreate_vers == 1
        const hid_t particle_index_dataset = H5Dcreate(h5_file, "id", HDF5_PARTICLE_INDEX_TYPE, index_dataspace, H5P_DEFAULT);
#else
        const hid_t particle_index_dataset = H5Dcreate(h5_file, "id", HDF5_PARTICLE_INDEX_TYPE, index_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

        // Select the local hyperslab from the dataspace
        H5Sselect_hyperslab(index_dataspace, H5S_SELECT_SET, offset, NULL, local_dataset_size, NULL);

        // Write index data to the HDF5 file
        H5Dwrite(particle_index_dataset, HDF5_PARTICLE_INDEX_TYPE, local_index_dataspace, index_dataspace, write_properties, &index_data[0]);

        H5Sclose(local_index_dataspace);
        H5Sclose(index_dataspace);
        H5Dclose(particle_index_dataset);

        // Write the property data
        for (unsigned int property = 0; property != property_component_list.size(); ++property)
          {
            const unsigned int data_components = (property_component_list[property].second != 2
                                                  ?
                                                  property_component_list[property].second
                                                  :
                                                  3);

            global_dataset_size[1] = data_components;
            local_dataset_size[1] = data_components;

            const hid_t property_dataspace = H5Screate_simple(2, global_dataset_size, NULL);
            const hid_t local_property_dataspace = H5Screate_simple(2, local_dataset_size, NULL);
#if H5Dcreate_vers == 1
            const hid_t particle_property_dataset = H5Dcreate(h5_file, property_component_list[property].first.c_str(), H5T_NATIVE_DOUBLE, property_dataspace, H5P_DEFAULT);
#else
            const hid_t particle_property_dataset = H5Dcreate(h5_file, property_component_list[property].first.c_str(), H5T_NATIVE_DOUBLE, property_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

            // Select the local hyperslab from the dataspace
            H5Sselect_hyperslab(property_dataspace, H5S_SELECT_SET, offset, NULL, local_dataset_size, NULL);

            // Write index data to the HDF5 file
            H5Dwrite(particle_property_dataset, H5T_NATIVE_DOUBLE, local_property_dataspace, property_dataspace, write_properties, &property_data[property][0]);

            H5Sclose(local_property_dataspace);
            H5Sclose(property_dataspace);
            H5Dclose(particle_property_dataset);
          }

        H5Pclose(write_properties);
        H5Fclose(h5_file);

        // Record and output XDMF info on root process
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            const std::string local_h5_filename = output_file_prefix + ".h5";
            XDMFEntry   entry(local_h5_filename, current_time, n_global_particles, 0, 3);
            DataOut<dim> data_out;
            const std::string xdmf_filename = (this->get_output_directory() + "particle.xdmf");

            entry.add_attribute("id", 1);

            for (unsigned int property = 0; property < property_component_list.size(); ++property)
              {
                const unsigned int data_components = (property_component_list[property].second != 2
                                                      ?
                                                      property_component_list[property].second
                                                      :
                                                      3);

                entry.add_attribute(property_component_list[property].first, data_components);
              }

            xdmf_entries.push_back(entry);

            data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), this->get_mpi_communicator());
          }

        file_index++;

        return output_path_prefix;
#else
        (void) property_component_list;
        (void) particles;
        (void) current_time;
        return "";
#endif
      }


      template <int dim>
      template <class Archive>
      void HDF5Output<dim>::serialize (Archive &ar, const unsigned int)
      {
        // invoke serialization of the base class
        ar &file_index
        &xdmf_entries
        ;
      }

      template <int dim>
      void
      HDF5Output<dim>::save (std::ostringstream &os) const
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      template <int dim>
      void
      HDF5Output<dim>::load (std::istringstream &is)
      {
        aspect::iarchive ia (is);
        ia >> (*this);
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

