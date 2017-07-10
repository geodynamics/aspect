/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <aspect/particle/output/hdf5.h>
#include <aspect/utilities.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/base/utilities.h>

#ifdef DEAL_II_WITH_HDF5
#include <hdf5.h>
#endif

namespace aspect
{
  namespace Particle
  {
    namespace Output
    {

      // Define the hdf5 type of the particle index output
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
      {}



      template <int dim>
      void HDF5Output<dim>::initialize ()
      {
#ifndef DEAL_II_WITH_HDF5
        AssertThrow (false,
                     ExcMessage ("deal.ii was not compiled with HDF5 support, "
                                 "so HDF5 output is not possible. Please "
                                 "recompile deal.ii with HDF5 support turned on "
                                 "or select a different particle output format."));
#endif

        aspect::Utilities::create_directory (this->get_output_directory() + "particles/",
                                             this->get_mpi_communicator(),
                                             true);
      }



      template <int dim>
      std::string
      HDF5Output<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                            const Property::ParticlePropertyInformation &property_information,
                                            const double current_time)
      {
#ifdef DEAL_II_WITH_HDF5
        // Create the filename
        const std::string output_file_prefix = "particles-" + Utilities::int_to_string (file_index, 5);
        const std::string output_path_prefix =
          this->get_output_directory()
          + "particles/"
          + output_file_prefix;
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
        MPI_Scan(&n_local_particles, &local_particle_index_offset, 1, ASPECT_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());

        hsize_t offset[2];
        offset[0] = local_particle_index_offset - n_local_particles;
        offset[1] = 0;

        // Prepare the output data
        std::vector<double> position_data (3 * n_local_particles,0.0);
        std::vector<types::particle_index> index_data (n_local_particles);
        std::vector<std::vector<double> > property_data;

        // Resize the property vectors
        for (unsigned int property = 0; property < property_information.n_fields(); ++property)
          {
            const unsigned int n_components = property_information.get_components_by_field_index(property);

            if (n_components == dim)
              property_data.push_back(std::vector<double> (3 * n_local_particles,0.0));
            else
              for (unsigned int component = 0; component < n_components; ++component)
                property_data.push_back(std::vector<double> (1 * n_local_particles));
          }

        // Write into the output vectors
        typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator it = particles.begin();
        for (unsigned int i = 0; it != particles.end(); ++i, ++it)
          {
            for (unsigned int d = 0; d < dim; ++d)
              position_data[i*3+d] = it->second.get_location()(d);

            index_data[i] = it->second.get_id();

            const ArrayView<const double> properties = it->second.get_properties();

            unsigned int particle_property_index = 0;

            unsigned int output_field_index = 0;
            for (unsigned int property = 0; property < property_information.n_fields(); ++property)
              {
                const unsigned int n_components = property_information.get_components_by_field_index(property);
                if (n_components == dim)
                  {
                    for (unsigned int component = 0; component < n_components; ++component,++particle_property_index)
                      property_data[output_field_index][i * 3 + component] = properties[particle_property_index];

                    ++output_field_index;
                  }
                else
                  for (unsigned int component = 0; component < n_components; ++component,++particle_property_index,++output_field_index)
                    property_data[output_field_index][i] = properties[particle_property_index];
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

        global_dataset_size[1] = 1;
        local_dataset_size[1] = 1;
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
        unsigned int output_index = 0;
        for (unsigned int property = 0; property != property_information.n_fields(); ++property)
          {
            const unsigned int n_components = property_information.get_components_by_field_index(property);
            const unsigned int data_components = ((n_components == dim)
                                                  ?
                                                  3
                                                  :
                                                  1);
            const unsigned int data_fields = ((n_components == dim) || (n_components == 1)) ?
                                             1
                                             :
                                             n_components;

            global_dataset_size[1] = data_components;
            local_dataset_size[1] = data_components;

            for (unsigned int i = 0; i < data_fields; ++i)
              {
                std::string field_name = property_information.get_field_name_by_index(property);
                if (data_fields > 1)
                  field_name.append('_' + Utilities::to_string(i));

                const hid_t property_dataspace = H5Screate_simple(2, global_dataset_size, NULL);
                const hid_t local_property_dataspace = H5Screate_simple(2, local_dataset_size, NULL);
#if H5Dcreate_vers == 1
                const hid_t particle_property_dataset = H5Dcreate(h5_file, field_name.c_str(), H5T_NATIVE_DOUBLE, property_dataspace, H5P_DEFAULT);
#else
                const hid_t particle_property_dataset = H5Dcreate(h5_file, field_name.c_str(), H5T_NATIVE_DOUBLE, property_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

                // Select the local hyperslab from the dataspace
                H5Sselect_hyperslab(property_dataspace, H5S_SELECT_SET, offset, NULL, local_dataset_size, NULL);

                // Write index data to the HDF5 file
                H5Dwrite(particle_property_dataset, H5T_NATIVE_DOUBLE, local_property_dataspace, property_dataspace, write_properties, &property_data[output_index][0]);

                H5Sclose(local_property_dataspace);
                H5Sclose(property_dataspace);
                H5Dclose(particle_property_dataset);
                ++output_index;
              }
          }

        H5Pclose(write_properties);
        H5Fclose(h5_file);

        // Record and output XDMF info on root process
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            const std::string local_h5_filename =
              "particles/"
              + output_file_prefix
              + ".h5";
            XDMFEntry   entry(local_h5_filename, current_time, n_global_particles, 0, 3);
            DataOut<dim> data_out;
            const std::string xdmf_filename = (this->get_output_directory() + "particles.xdmf");

            entry.add_attribute("id", 1);

            for (unsigned int property = 0; property < property_information.n_fields(); ++property)
              {
                const unsigned int n_components = property_information.get_components_by_field_index(property);

                if (n_components == 1)
                  entry.add_attribute(property_information.get_field_name_by_index(property), 1);
                else if (n_components == dim)
                  entry.add_attribute(property_information.get_field_name_by_index(property), 3);
                else
                  for (unsigned int component = 0; component < n_components; ++component)
                    entry.add_attribute(property_information.get_field_name_by_index(property) + "_" + Utilities::to_string(component), 1);
              }

            xdmf_entries.push_back(entry);

            data_out.write_xdmf_file(xdmf_entries, xdmf_filename.c_str(), this->get_mpi_communicator());
          }

        file_index++;

        return output_path_prefix;
#else
        (void) property_information;
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
        ar &static_cast<Interface<dim> &>(*this);

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

