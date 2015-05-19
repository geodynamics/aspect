/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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
          /**
           * Constructor.
           *
           * @param[in] The directory into which output files shall be placed.
           * @param[in] The MPI communicator that describes this simulation.
           */
        template <int dim>
        ASCIIOutput<dim>::ASCIIOutput()
                      :
                      Interface<dim> ()
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
        template <int dim>
          std::string
          ASCIIOutput<dim>::output_particle_data(const std::multimap<LevelInd, BaseParticle<dim> > &particles,
                               std::vector<MPIDataInfo> &data_info,
                               const double &current_time)
          {
            typename std::multimap<LevelInd, BaseParticle<dim> >::const_iterator  it;
            unsigned int                            i;
            std::string                             output_file_prefix, output_path_prefix, full_filename;
            std::vector<MPIDataInfo>::iterator      dit;

            output_file_prefix = "particle-" + Utilities::int_to_string (this->file_index, 5);
            output_path_prefix = this->output_dir + output_file_prefix;
            full_filename = output_path_prefix + "." + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->communicator), 4) + ".txt";
            std::ofstream output (full_filename.c_str());
            if (!output)
              std::cout << "ERROR: proc " << Utilities::MPI::this_mpi_process(this->communicator) << " could not create " << full_filename << std::endl;

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
                                               "")
    }
  }
}

