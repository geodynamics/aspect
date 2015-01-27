/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_output_h
#define __aspect__particle_output_h

#include <deal.II/base/mpi.h>
#include <aspect/particle/particle.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      /**
       * Abstract base class used for classes that generate particle output
       */
      template <int dim, class T>
      class Interface
      {
        protected:
          /**
           * Path to directory in which to put particle output files
           */
          const std::string     output_dir;

          /**
           * MPI communicator to be used for output synchronization
           */
          MPI_Comm        communicator;

//TODO: This needs to be serialized
          /**
           * Internal index of file output number, must be incremented by
           * derived classes when they create a new file.
           */
          unsigned int    file_index;

        public:
          /**
           * Constructor.
           *
           * @param[in] output_directory The directory into which output files
           * shall be placed.
           * @param[in] communicator The MPI communicator that describes this
           * simulation.
           */
          Interface(const std::string &output_directory,
                    const MPI_Comm     communicator)
            :
            output_dir (output_directory),
            communicator (communicator),
            file_index (0)
          {}

          /**
           * Destructor. Made virtual so that derived classes can be created
           * and destroyed through pointers to the base class.
           */
          virtual ~Interface () {}

          /**
           * Write data about the particles specified in the first argument to
           * a file. If possible, encode the current simulation time into this
           * file using the data provided in the second argument.
           *
           * @param [in] particles The set of particles to generate a
           * graphical representation for
           * @param [in] current_time Current time of the simulation, given as
           * either years or seconds, as selected in the input file. In other
           * words, output writers do not need to know the units in which time
           * is described.
           * @return The name of the file that was written, or any other
           * information that describes what output was produced if for
           * example multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<LevelInd, T> &particles,
                               const double &current_time) = 0;

          /**
           * Read or write the data of this object for serialization
           */
          template <class Archive>
          void serialize(Archive &ar, const unsigned int version)
          {
            ar &file_index
            ;
          }
      };


      /**
       * Create an output object.
       *
       * @param[in] data_format_name Name of the format in which the created
       * output writer should produce its files
       * @param[in] output_directory Directory into which to put the data
       * files
       * @param[in] communicator MPI communicator object that describes this
       * simulation
       * @return a pointer to an output object, needs to be deleted by the
       * caller.
       */
      template <int dim, class T>
      Interface<dim, T> *
      create_output_object (const std::string &data_format_name,
                            const std::string &output_directory,
                            const MPI_Comm     communicator);


      /**
       * Return a list of names (separated by '|') of possible writers of
       * graphical output formats for particle data.
       */
      std::string
      output_object_names ();

    }
  }
}

#endif
