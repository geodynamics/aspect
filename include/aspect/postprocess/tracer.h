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

#ifndef __aspect__postprocess_tracer_h
#define __aspect__postprocess_tracer_h

#include <aspect/postprocess/interface.h>
#include <aspect/particle/world.h>
#include <aspect/particle/output.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class PassiveTracers : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      private:
        // The world holding the particles
        Particle::World<dim, Particle::BaseParticle<dim> >              _world;

        // The integrator to use in moving the particles
        Particle::Integrator<dim, Particle::BaseParticle<dim> >         *_integrator;

        // Abstract output object
        Particle::Output<dim, Particle::BaseParticle<dim> >             *_output;

        // Whether this set has been initialized yet or not
        bool                            _initialized;

        // Number of initial particles to create
        // Use a double rather than int since doubles can represent up to 2^52
        double                          _num_initial_tracers;

        // Interval between output (in years if appropriate
        // simulation parameter is set, otherwise seconds)
        double                          _data_output_interval;

        // Output format for particle data
        std::string                     _data_output_format;

        // Integration scheme to move particles
        std::string                     _integration_scheme;

        // Records time for next output to occur
        double                          _next_data_output_time;

      public:
        PassiveTracers(void) : _initialized(false), _next_data_output_time(std::numeric_limits<double>::quiet_NaN()) {};

        virtual std::pair<std::string,std::string> execute (TableHandler &statistics);

        void set_next_data_output_time (const double current_time);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
  }
}

#endif
