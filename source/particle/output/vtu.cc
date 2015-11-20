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

#include <aspect/particle/output/vtu.h>

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
      std::string
      VTUOutput<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                           const std::vector<std::pair<std::string, unsigned int> > &property_component_list,
                                           const double current_time)
      {
        const std::string output_file_prefix = "particles-" + Utilities::int_to_string (file_index, 5);
        const std::string output_path_prefix = this->get_output_directory() + output_file_prefix;

        const std::string filename = (output_file_prefix +
                                      "." +
                                      Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()), 4) +
                                      ".vtu");
        const std::string full_filename = (this->get_output_directory() + filename);

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
        std::vector<std::pair<std::string,unsigned int> >::const_iterator property = property_component_list.begin();
        unsigned int data_offset = 0;

        for (; property!=property_component_list.end(); ++property)
          {

            output << "        <DataArray type=\"Float64\" Name=\"" << property->first << "\" NumberOfComponents=\"" << (property->second == 2 ? 3 : property->second) << "\" Format=\"ascii\">\n";
            for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator
                 it=particles.begin(); it!=particles.end(); ++it)
              {
                const std::vector<double> particle_data = it->second.get_properties();

                output << "          ";

                for (unsigned int d=0; d < property->second; ++d)
                  {
                    output << particle_data[data_offset+d];

                    if (d+1 < property->second)
                      output << ' ';
                  }

                if (property->second == 2)
                  output << " 0";
                output << "\n";
              }
            data_offset += property->second;
            output << "        </DataArray>\n";
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
            const std::string full_pvtu_filename = (this->get_output_directory() + pvtu_filename);

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

            for (property=property_component_list.begin(); property!=property_component_list.end(); ++property)
              {
                pvtu_output << "      <PDataArray type=\"Float64\" Name=\"" << property->first << "\" NumberOfComponents=\"" << (property->second == 2 ? 3 : property->second) << "\" format=\"ascii\"/>\n";
              }
            pvtu_output << "    </PPointData>\n";
            for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++i)
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
            const std::string pvd_master_filename = (this->get_output_directory() + "particles.pvd");
            std::ofstream pvd_master (pvd_master_filename.c_str());
            DataOut<dim>().write_pvd_record (pvd_master, times_and_pvtu_file_names);

//TODO: write a global .visit record. this needs a variant of the write_visit_record
// function
            /*
                          const std::string visit_master_filename = (this->get_output_directory() + "particles.visit");
                          std::ofstream visit_master (visit_master_filename.c_str());
                          DataOut<dim>().write_visit_record (visit_master, vtu_file_names);
            */
          }
        file_index++;

        return output_path_prefix;
      }

      template <int dim>
      template <class Archive>
      void VTUOutput<dim>::serialize (Archive &ar, const unsigned int)
      {
        // invoke serialization of the base class
        ar &file_index
        & times_and_pvtu_file_names
        & vtu_file_names
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


