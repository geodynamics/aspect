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

#ifndef __aspect__postprocess_tracer_h
#define __aspect__postprocess_tracer_h

#include <aspect/postprocess/interface.h>
#include <aspect/particle/world.h>
#include <aspect/particle/output/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/particle/particle.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A Postprocessor that creates tracer particles, which follow the
     * velocity field of the simulation. The particles can be generated
     * and propagated in various ways and they can carry a number of
     * constant or time-varying properties. The postprocessor can write
     * output positions and properties of all tracers at chosen intervals,
     * although this is not mandatory. It also allows other parts of the
     * code to query the tracers for information.
     */
    template <int dim>
    class Tracers : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Tracers();

        /**
         * Destructor.
         */
        ~Tracers();

        /**
         *
         */
        void
        initialize();

        /**
         * Returns a const reference to the particle world, in case anyone
         * wants to query something about tracers.
         */
        const Particle::World<dim> &
        get_particle_world() const;

        /**
         * Returns a reference to the particle world, in case anyone wants to
         * change something within the particle world. Use with care, usually
         * you want to only let the functions within the particle subsystem
         * change member variables of the particle world.
         */
        Particle::World<dim> &
        get_particle_world();

        /**
         * Execute this postprocessor. Derived classes will implement this
         * function to do whatever they want to do to evaluate the solution at
         * the current time step.
         *
         * @param[in,out] statistics An object that contains statistics that
         * are collected throughout the simulation and that will be written to
         * an output file at the end of each time step. Postprocessors may
         * deposit data in these tables for later visualization or further
         * processing.
         *
         * @return A pair of strings that will be printed to the screen after
         * running the postprocessor in two columns; typically the first
         * column contains a description of what the data is and the second
         * contains a numerical value of this data. If there is nothing to
         * print, simply return two empty strings.
         */
        virtual
        std::pair<std::string,std::string> execute (TableHandler &statistics);

        /**
         * Save the state of this object.
         */
        virtual
        void save (std::map<std::string, std::string> &status_strings) const;

        /**
         * Restore the state of the object.
         */
        virtual
        void load (const std::map<std::string, std::string> &status_strings);

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * The world holding the particles
         */
        Particle::World<dim> world;

        /**
         * Pointer to an output object
         */
        std_cxx11::unique_ptr<Particle::Output::Interface<dim> >       output;

        /**
         * Whether particles have been created and initialized yet or not.
         * We can not simply use the number of particles in the domain to
         * check this, because all particles could leave the domain through an
         * open boundary. This variable is initialized to false and set to
         * true in the first call of execute().
         */
        bool                            initialized;

        /**
         * Interval between output (in years if appropriate simulation
         * parameter is set, otherwise seconds)
         */
        double                          data_output_interval;

        /**
         * Records time for next output to occur
         */
        double                          next_data_output_time;

        /**
         * Compute the next time output should be generated assuming that
         * output was generated at current_time, and set the result in the
         * next_data_output_time variable.
         */
        void set_next_data_output_time (const double current_time);
    };
  }
}

#endif
