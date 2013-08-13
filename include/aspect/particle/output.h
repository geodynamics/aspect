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
/*  $Id$  */

#ifndef __aspect__particle_output_h
#define __aspect__particle_output_h

#include <deal.II/base/mpi.h>
#include <aspect/particle/particle.h>


namespace aspect
{
  namespace Particle
  {
    /**
     *  Abstract base class used for classes that generate particle output
     */
    template <int dim, class T>
    class Output
    {
      protected:
        /**
         *  MPI communicator to be used for output synchronization
         */
        MPI_Comm        communicator;

        /**
         *  Internal index of file output number, must be incremented
         *  by derived classes when they create a new file.
         */
        unsigned int    file_index;


        // Path to directory in which to put particle output files
        std::string     output_dir;

      public:
        /**
         * Constructor.
         */
        Output()
      :
        file_index (0)
        {}

        unsigned int self_rank()
        {
          return Utilities::MPI::this_mpi_process(communicator);
        };

        unsigned int world_size()
        {
          return Utilities::MPI::n_mpi_processes(communicator);
        };

        void set_mpi_comm(MPI_Comm new_comm_world)
        {
          communicator = new_comm_world;
        };

        void set_output_directory(const std::string &new_out_dir)
        {
          output_dir = new_out_dir;
        };

        virtual std::string output_particle_data(const std::multimap<LevelInd, T> &particles,
                                                 const double &current_time) = 0;
    };


    /**
     * Create an output object given the specified name.
     */
    template <int dim, class T>
    Output<dim, T> *
    create_output_object (const std::string &data_format_name);


    /**
     * Return a list of names (separated by '|') of possible writers of graphical
     * output formats for particle data.
     */
    std::string
    output_object_names ();
  }
}

#endif
