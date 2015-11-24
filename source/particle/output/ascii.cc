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

#include <aspect/particle/output/ascii.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      template <int dim>
      ASCIIOutput<dim>::ASCIIOutput()
        :
        file_index(0)
      {}

      template <int dim>
      std::string
      ASCIIOutput<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                             const std::vector<std::pair<std::string, unsigned int> > &property_component_list,
                                             const double /*time*/)
      {
        const std::string output_file_prefix = "particle-" + Utilities::int_to_string (file_index, 5);
        const std::string output_path_prefix = this->get_output_directory() + output_file_prefix;
        const std::string full_filename = output_path_prefix + "." + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()), 4) + ".txt";

        std::ofstream output (full_filename.c_str());

        AssertThrow (output,
                     ExcMessage (std::string("Could not open ascii particle output file <"
                                             +
                                             full_filename
                                             +
                                             ">.")));

        // Print the header line
        output << "# ";
        for (unsigned int i = 0; i < dim; ++i)
          output << "position[" << i << "] ";

        output << "id ";

        std::vector<std::pair<std::string,unsigned int> >::const_iterator property = property_component_list.begin();
        for (; property!=property_component_list.end(); ++property)
          {
            // If it is a 1D element, print just the name, otherwise use []
            if (property->second == 1)
              output << property->first << ' ';
            else
              for (unsigned int component_index=0; component_index<property->second; ++component_index)
                output << property->first << "[" << component_index << "] ";
          }
        output << "\n";

        // And print the data for each particle
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator it=particles.begin(); it!=particles.end(); ++it)
          {
            const std::vector<double> properties = it->second.get_properties();

            output << it->second.get_location();
            output << ' ' << it->second.get_id();
            for (unsigned int i = 0; i < properties.size(); ++i)
              output << ' ' << properties[i];

            output << "\n";
          }

        file_index++;

        return output_path_prefix;
      }

      template <int dim>
      template <class Archive>
      void ASCIIOutput<dim>::serialize (Archive &ar, const unsigned int)
      {
        // invoke serialization of the base class
        ar &file_index
        ;
      }

      template <int dim>
      void
      ASCIIOutput<dim>::save (std::ostringstream &os) const
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      template <int dim>
      void
      ASCIIOutput<dim>::load (std::istringstream &is)
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
      ASPECT_REGISTER_PARTICLE_OUTPUT(ASCIIOutput,
                                      "ascii",
                                      "This particle output plugin writes particle "
                                      "positions and properties into space separated "
                                      "ascii files.")
    }
  }
}

