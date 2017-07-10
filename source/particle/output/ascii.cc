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

#include <aspect/particle/output/ascii.h>
#include <aspect/utilities.h>


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
      void ASCIIOutput<dim>::initialize ()
      {
        aspect::Utilities::create_directory (this->get_output_directory() + "particles/",
                                             this->get_mpi_communicator(),
                                             true);
      }

      template <int dim>
      std::string
      ASCIIOutput<dim>::output_particle_data(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                             const Property::ParticlePropertyInformation &property_information,
                                             const double /*time*/)
      {
        const std::string output_file_prefix =
          this->get_output_directory()
          + "particles/"
          + "particle-"
          + Utilities::int_to_string (file_index, 5);
        const std::string full_filename =
          output_file_prefix
          + "."
          + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()), 4)
          + ".txt";

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

        output << "id";

        for (unsigned int field_index = 0; field_index < property_information.n_fields(); ++field_index)
          {
            const unsigned n_components = property_information.get_components_by_field_index(field_index);
            const std::string field_name = property_information.get_field_name_by_index(field_index);
            // If it is a 1D element, print just the name, otherwise use []
            if (n_components == 1)
              output << ' ' << field_name;
            else
              for (unsigned int component_index=0; component_index<n_components; ++component_index)
                output << ' ' << field_name << "[" << component_index << "]";
          }
        output << "\n";

        // And print the data for each particle
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator it=particles.begin(); it!=particles.end(); ++it)
          {

            output << it->second.get_location();
            output << ' ' << it->second.get_id();

            if (property_information.n_fields() > 0)
              {
                const ArrayView<const double> properties = it->second.get_properties();

                for (unsigned int i = 0; i < properties.size(); ++i)
                  output << ' ' << properties[i];
              }

            output << "\n";
          }

        file_index++;

        return output_file_prefix;
      }

      template <int dim>
      template <class Archive>
      void ASCIIOutput<dim>::serialize (Archive &ar, const unsigned int)
      {
        // invoke serialization of the base class
        ar &static_cast<Interface<dim> &>(*this);

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

