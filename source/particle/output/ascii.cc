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
        Interface<dim> ()
      {}

      template <int dim>
      std::string
      ASCIIOutput<dim>::output_particle_data(const std::multimap<LevelInd, Particle<dim> > &particles,
                                             const std::vector<std::string> &names,
                                             const std::vector<unsigned int> &lengths,
                                             const double &)
      {
        const std::string output_file_prefix = "particle-" + Utilities::int_to_string (this->file_index, 5);
        const std::string output_path_prefix = this->output_dir + output_file_prefix;
        const std::string full_filename = output_path_prefix + "." + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->communicator), 4) + ".txt";
        std::ofstream output (full_filename.c_str());
        if (!output)
          std::cout << "ERROR: proc " << Utilities::MPI::this_mpi_process(this->communicator) << " could not create " << full_filename << std::endl;

        // Print the header line
        output << "# ";
        for (unsigned int i = 0; i < dim; ++i)
          output << "position[" << i << "] ";

        output << "id ";

        std::vector<std::string>::const_iterator  name = names.begin();
        std::vector<unsigned int>::const_iterator length = lengths.begin();
        for (; name!=names.end(),length!=lengths.end(); ++name,++length)
          {
            // If it's a 1D element, print just the name, otherwise use []
            if (*length == 1)
              {
                output << *name << " ";
              }
            else
              {
                for (unsigned int i=0; i<*length; ++i) output << *name << "[" << i << "] ";
              }
          }
        output << "\n";

        // And print the data for each particle
        for (typename std::multimap<LevelInd, Particle<dim> >::const_iterator it=particles.begin(); it!=particles.end(); ++it)
          {
            std::vector<double>  particle_data;
            it->second.write_data(particle_data);
            for (unsigned int i = 0; i < particle_data.size(); ++i)
              {
                output << particle_data[i] << " ";
              }
            output << "\n";
          }

        output.close();

        this->file_index++;

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
      ASPECT_REGISTER_PARTICLE_OUTPUT(ASCIIOutput,
                                      "ascii",
                                      "This particle output plugin writes particle "
                                      "positions and properties into space separated "
                                      "ascii files.")
    }
  }
}

