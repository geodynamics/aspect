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

#ifndef __aspect__particle_output_interface_h
#define __aspect__particle_output_interface_h

#include <aspect/particle/base_particle.h>
#include <aspect/particle/definitions.h>
#include <aspect/plugins.h>

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/mpi.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      /**
       * Abstract base class used for classes that generate particle output
       */
      template <int dim>
      class Interface
      {
        protected:
          /**
           * Path to directory in which to put particle output files
           */
          std::string     output_dir;

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
           * @param[in] The directory into which output files shall be placed.
           * @param[in] The MPI communicator that describes this simulation.
           */
          Interface() {}

          /**
           * Destructor. Made virtual so that derived classes can be created
           * and destroyed through pointers to the base class.
           */
          virtual ~Interface () {}


          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run and after the
           * SimulatorAccess (if applicable) is initialized.
           *
           * @param[in] output_directory The directory into which output files
           * shall be placed.
           * @param[in] communicator The MPI communicator that describes this
           * simulation.
           */
          virtual
          void
          initialize (std::string output_directory,
                      MPI_Comm communicator_);

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
           * @return The name of the file that was written, or
           * any other information that describes what output was produced if
           * for example multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<LevelInd, BaseParticle<dim> > &particles,
                               std::vector<MPIDataInfo> &data_info,
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

          /**
           * Declare the parameters this class takes through input files. The
           * default implementation of this function does not describe any
           * parameters. Consequently, derived classes do not have to overload
           * this function if they do not take any runtime parameters.
           */
          static
          void
          declare_parameters (ParameterHandler &);

          /**
           * Read the parameters this class declares from the parameter file.
           * The default implementation of this function does not read any
           * parameters. Consequently, derived classes do not have to overload
           * this function if they do not take any runtime parameters.
           */
          virtual
          void
          parse_parameters (ParameterHandler &);
      };


      /**
       * Register a particle output so that it can be selected from
       * the parameter file.
       *
       * @param name A string that identifies the particle output
       * @param description A text description of what this output does and that
       * will be listed in the documentation of the parameter file.
       * @param declare_parameters_function A pointer to a function that can be
       * used to declare the parameters that this particle output wants to read
       * from input files.
       * @param factory_function A pointer to a function that can create an
       * object of this particle output.
       *
       * @ingroup ParticleOutputs
       */
      template <int dim>
      void
      register_particle_output (const std::string &name,
                                     const std::string &description,
                                     void (*declare_parameters_function) (ParameterHandler &),
                                     Interface<dim> *(*factory_function) ());

      /**
       * A function that given the name of a model returns a pointer to an
       * object that describes it. Ownership of the pointer is transferred to
       * the caller.
       *
       * The model object returned is not yet initialized and has not
       * read its runtime parameters yet.
       *
       * @ingroup ParticleOutputs
       */
      template <int dim>
      Interface<dim> *
      create_particle_output (ParameterHandler &prm);

      /**
       * Declare the runtime parameters of the registered particle outputs.
       *
       * @ingroup ParticleOutputs
       */
      template <int dim>
      void
      declare_parameters (ParameterHandler &prm);

/**
 * Given a class name, a name, and a description for the parameter file
 * for a particle output, register it with the functions that
 * can declare their parameters and create these objects.
 *
 * @ingroup ParticleOutputs
 */
#define ASPECT_REGISTER_PARTICLE_OUTPUT(classname, name, description) \
template class classname<2>; \
template class classname<3>; \
namespace ASPECT_REGISTER_PARTICLE_OUTPUT_ ## classname \
{ \
aspect::internal::Plugins::RegisterHelper<aspect::Particle::Output::Interface<2>,classname<2> > \
dummy_ ## classname ## _2d (&aspect::Particle::Output::register_particle_output<2>, \
                            name, description); \
aspect::internal::Plugins::RegisterHelper<aspect::Particle::Output::Interface<3>,classname<3> > \
dummy_ ## classname ## _3d (&aspect::Particle::Output::register_particle_output<3>, \
                            name, description); \
}
}
}
}

#endif
