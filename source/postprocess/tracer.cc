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

#include <aspect/global.h>
#include <aspect/postprocess/tracer.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    PassiveTracers<dim>::PassiveTracers ()
      :
      integrator(NULL),
      output(NULL),
      generator(NULL),
      property_manager(),
      initialized(false),
      next_data_output_time(std::numeric_limits<double>::quiet_NaN())
    {}


    template <int dim>
    PassiveTracers<dim>::~PassiveTracers ()
    {
      if (integrator) delete integrator;
      if (output) delete output;
      if (generator) delete generator;
    }

    template <int dim>
    void
    PassiveTracers<dim>::initialize()
    {

    }

    template <int dim>
    const Particle::World<dim> &
    PassiveTracers<dim>::get_particle_world() const
    {
      return world;
    }

    template <int dim>
    std::pair<std::string,std::string>
    PassiveTracers<dim>::execute (TableHandler &statistics)
    {
      if (!initialized)
        {
          next_data_output_time = this->get_time();

          // Set up the particle world with the appropriate simulation objects
          world.set_integrator(integrator);
          world.set_manager(&property_manager);

          // And initialize the world
          world.init();

          next_data_output_time = this->get_time();

          // Add the specified number of particles
          generator->generate_particles(world);
          world.initialize_particles();

          initialized = true;
        }

      // Advance the particles in the world to the current time
      world.advance_timestep();

      const unsigned int num_particles = world.get_global_particle_count();
      statistics.add_value("Advected particles",num_particles);
      std::ostringstream result_string;
      result_string << num_particles;

      // If it's time to generate an output file, call the appropriate functions and reset the timer
      if (this->get_time() >= next_data_output_time)
        {
          if (property_manager.need_update() == Particle::Property::update_output_step)
            world.update_particles();

          set_next_data_output_time (this->get_time());

          std::vector<std::string> names;
          std::vector<unsigned int> length;

          property_manager.get_data_info(names,
                                         length);

          const std::string data_file_name = output->output_particle_data(world.get_particles(),
                                                                          names,
                                                                          length,
                                                                          (this->convert_output_to_years() ?
                                                                           this->get_time() / year_in_seconds :
                                                                           this->get_time()));

          result_string << ". Writing particle graphical output " + data_file_name;
        }

      return std::make_pair("Advected particles: ", result_string.str());
    }



    template <int dim>
    void
    PassiveTracers<dim>::set_next_data_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
      if (data_output_interval > 0)
        {
          // the current time is always in seconds, so we need to convert the output_interval to the same unit
          const double output_interval_in_s = (this->convert_output_to_years() ?
                                               (data_output_interval*year_in_seconds) :
                                               data_output_interval);

          // we need to compute the smallest integer that is bigger than current_time/my_output_interval,
          // even if it is a whole number already (otherwise we output twice in a row)
          next_data_output_time = (std::floor(current_time/output_interval_in_s)+1.0) * output_interval_in_s;
        }
    }

    template <int dim>
    void
    PassiveTracers<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Tracers");
        {
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step.\n\n"
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      Particle::Generator::declare_parameters<dim>(prm);
      Particle::Output::declare_parameters<dim>(prm);
      Particle::Integrator::declare_parameters<dim>(prm);
      Particle::Property::Manager<dim>::declare_parameters(prm);
    }


    template <int dim>
    void
    PassiveTracers<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Tracers");
        {
          data_output_interval = prm.get_double ("Time between data output");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      world.initialize(this->get_simulator());

      // Create a generator object using a random uniform distribution
      generator = Particle::Generator::create_particle_generator<dim>
                  (prm);
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(generator))
        sim->initialize (this->get_simulator());
      generator->parse_parameters(prm);

      // Create an output object depending on what the parameters specify
      output = Particle::Output::create_particle_output<dim>
               (prm);
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(output))
        sim->initialize (this->get_simulator());
      output->initialize(this->get_output_directory(),
                         this->get_mpi_communicator());

      // Create an integrator object depending on the specified parameter
      integrator = Particle::Integrator::create_particle_integrator<dim>
                   (prm);
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(integrator))
        sim->initialize (this->get_simulator());

      SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&property_manager);
      sim->initialize (this->get_simulator());
      property_manager.parse_parameters(prm);
      property_manager.initialize();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PassiveTracers,
                                  "tracers",
                                  "Postprocessor that propagates passive tracer particles based on the "
                                  "velocity field.")
  }
}
