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

#include <aspect/particle/output/vtu.h>
#include <aspect/utilities.h>

#include <deal.II/numerics/data_out.h>

namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      template <int dim>
      VTUOutput<dim>::VTUOutput()
        :
        file_index(0)
      {}

      template <int dim>
      void VTUOutput<dim>::initialize ()
      {
        aspect::Utilities::create_directory (this->get_output_directory() + "particles/",
                                             this->get_mpi_communicator(),
                                             true);
      }

      template <int dim>
      std::string
      VTUOutput<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                           const Property::ParticlePropertyInformation &property_information,
                                           const double current_time)
      {
        const std::string output_file_prefix = "particles-" + Utilities::int_to_string (file_index, 5);
        const std::string output_path_prefix = this->get_output_directory()
                                               + "particles/"
                                               + output_file_prefix;

        const std::string filename = (output_file_prefix +
                                      "." +
                                      Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()), 4) +
                                      ".vtu");
        const std::string full_filename = this->get_output_directory()
                                          + "particles/"
                                          + filename;

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
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator
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

        // Write data for each particle (id, velocity, etc)
        output << "      <PointData Scalars=\"scalars\">\n";

        output << "        <DataArray type=\"UInt64\" Name=\"id\" NumberOfComponents=\"1\" Format=\"ascii\">\n";
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator
             it=particles.begin(); it!=particles.end(); ++it)
          output << "          " << it->second.get_id() << "\n" ;

        output << "        </DataArray>\n";

        // Print the data associated with the particles
        unsigned int data_offset = 0;

        for (unsigned int field_index = 0; field_index < property_information.n_fields(); ++field_index)
          {
            const unsigned int n_components = property_information.get_components_by_field_index(field_index);

            // If this is a property with one component or as many components
            // as spatial dimensions, output it as one scalar / vector field.
            // Vector fields are padded with zeroes to 3 dimensions (vtk limitation).
            if ((n_components == 1) || (n_components == dim))
              {
                output << "        <DataArray type=\"Float64\" Name=\""
                       << property_information.get_field_name_by_index(field_index)
                       << "\" NumberOfComponents=\"" << (n_components == 1 ? 1 : 3)
                       << "\" Format=\"ascii\">\n";
                for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator
                     it=particles.begin(); it!=particles.end(); ++it)
                  {
                    const ArrayView<const double> particle_data = it->second.get_properties();

                    output << "         ";
                    for (unsigned int d=0; d < n_components; ++d)
                      output << ' ' << particle_data[data_offset+d];

                    if (n_components == 2)
                      output << " 0";
                    output << "\n";
                  }
                data_offset += n_components;
                output << "        </DataArray>\n";
              }
            // Otherwise create n_components scalar fields
            else
              {
                for (unsigned int d=0; d<n_components; ++d)
                  {
                    output << "        <DataArray type=\"Float64\" Name=\""
                           << property_information.get_field_name_by_index(field_index) << '_' << d
                           << "\" NumberOfComponents=\"1"
                           << "\" Format=\"ascii\">\n";
                    for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator
                         it=particles.begin(); it!=particles.end(); ++it)
                      {
                        const ArrayView<const double> particle_data = it->second.get_properties();

                        output << particle_data[data_offset+d] << "\n";
                      }
                    output << "        </DataArray>\n";
                  }
                data_offset += n_components;

              }
          }
        output << "      </PointData>\n";

        output << "    </Piece>\n";
        output << "  </UnstructuredGrid>\n";
        output << "</VTKFile>\n";

        output.close();


        // Write the parallel pvtu and pvd files on the root process
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            const std::string pvtu_filename = (output_file_prefix + ".pvtu");
            const std::string full_pvtu_filename = this->get_output_directory()
                                                   + "particles/"
                                                   + pvtu_filename;

            std::ofstream pvtu_output (full_pvtu_filename.c_str());
            AssertThrow (pvtu_output, ExcIO());

            pvtu_output << "<?xml version=\"1.0\"?>\n";
            pvtu_output << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
            pvtu_output << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
            pvtu_output << "    <PPoints>\n";
            pvtu_output << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
            pvtu_output << "    </PPoints>\n";
            pvtu_output << "    <PPointData Scalars=\"scalars\">\n";
            pvtu_output << "      <PDataArray type=\"UInt64\" Name=\"id\" NumberOfComponents=\"1\" Format=\"ascii\"/>\n";

            for (unsigned int field_index = 0; field_index < property_information.n_fields(); ++field_index)
              {
                const unsigned int n_components = property_information.get_components_by_field_index(field_index);
                if ((n_components == 1) || (n_components == dim))
                  pvtu_output << "      <PDataArray type=\"Float64\" Name=\"" << property_information.get_field_name_by_index(field_index)
                              << "\" NumberOfComponents=\"" << (n_components == 1 ? 1 : 3)
                              << "\" format=\"ascii\"/>\n";
                else
                  for (unsigned int component = 0; component < property_information.get_components_by_field_index(field_index); ++component)
                    pvtu_output << "      <PDataArray type=\"Float64\" Name=\"" << property_information.get_field_name_by_index(field_index)
                                << '_' << component
                                << "\" NumberOfComponents=\"1"
                                << "\" format=\"ascii\"/>\n";
              }
            pvtu_output << "    </PPointData>\n";
            for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++i)
              {
                pvtu_output << "    <Piece Source=\"" << output_file_prefix << "." << Utilities::int_to_string(i, 4) << ".vtu\"/>\n";
              }
            pvtu_output << "  </PUnstructuredGrid>\n";
            pvtu_output << "</VTKFile>\n";
            pvtu_output.close();

            // update the .pvd record for the entire simulation
            times_and_pvtu_file_names.push_back(std::make_pair(current_time,
                                                               "particles/"+pvtu_filename));
            const std::string pvd_master_filename = (this->get_output_directory() + "particles.pvd");
            std::ofstream pvd_master (pvd_master_filename.c_str());

            DataOutBase::write_pvd_record (pvd_master, times_and_pvtu_file_names);

            // same for the .visit record for the entire simulation. for this, we first
            // have to collect all files that together form this one time step
            std::vector<std::string> this_timestep_output_files;
            for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++i)
              this_timestep_output_files.push_back ("particles/" + output_file_prefix +
                                                    "." + Utilities::int_to_string(i, 4) + ".vtu");
            times_and_vtu_file_names.push_back (std::make_pair (current_time,
                                                                this_timestep_output_files));

            const std::string visit_master_filename = (this->get_output_directory() + "particles.visit");
            std::ofstream visit_master (visit_master_filename.c_str());
            DataOutBase::write_visit_record (visit_master, times_and_vtu_file_names);
          }
        file_index++;

        return output_path_prefix;
      }

      template <int dim>
      template <class Archive>
      void VTUOutput<dim>::serialize (Archive &ar, const unsigned int)
      {
        // invoke serialization of the base class
        ar &static_cast<Interface<dim> &>(*this);

        ar &file_index
        & times_and_pvtu_file_names
        & times_and_vtu_file_names
        ;
      }

      template <int dim>
      void
      VTUOutput<dim>::save (std::ostringstream &os) const
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      template <int dim>
      void
      VTUOutput<dim>::load (std::istringstream &is)
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
      ASPECT_REGISTER_PARTICLE_OUTPUT(VTUOutput,
                                      "vtu",
                                      "This particle output plugin writes particle "
                                      "positions and properties into vtu files.")
    }
  }
}


