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

#include <aspect/particle/output/null.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      // Blank class to avoid output of data, for example, if particles are used to represent
      // other physical phenomenon that influences the simulation and we don't care about their positions

      /**
       * Constructor.
       *
       * @param[in] The directory into which output files shall be placed.
       * @param[in] The MPI communicator that describes this simulation.
       */
      template <int dim>
      NullOutput<dim>::NullOutput()
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
      NullOutput<dim>::output_particle_data(const std::multimap<LevelInd, BaseParticle<dim> > &/*particles*/,
                                            const std::vector<std::string>  &/*data_names*/,
                                            const std::vector<unsigned int> &/*data_components*/,
                                            const double &/*current_time*/)
      {
        return "";
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
      ASPECT_REGISTER_PARTICLE_OUTPUT(NullOutput,
                                      "none",
                                      "Blank class to avoid output of data, for example, if "
                                      "particles are used to represent other physical phenomenon "
                                      "that influences the simulation and we don't care about "
                                      "their positions.")
    }
  }
}


