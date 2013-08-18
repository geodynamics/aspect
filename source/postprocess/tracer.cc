/*
 Copyright (C) 2011, 2012, 2013 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/postprocess/tracer.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    PassiveTracers<dim>::PassiveTracers ()
      :
      initialized(false),
      next_data_output_time(std::numeric_limits<double>::quiet_NaN())
    {}



    template <int dim>
    std::pair<std::string,std::string>
    PassiveTracers<dim>::execute (TableHandler &statistics)
    {
      bool            output_data = false;

      if (!initialized)
        {
          // Create an output object depending on what the parameters specify
          output = Particle::Output::create_output_object<dim,Particle::BaseParticle<dim> >
                   (data_output_format,
                    this->get_output_directory(),
                    this->get_mpi_communicator());

          //TODO: Use factory methods here
          // Create an integrator object depending on the specified parameter
          if (integration_scheme == "euler")
            integrator = new Particle::EulerIntegrator<dim, Particle::BaseParticle<dim> >;
          else if (integration_scheme == "rk2")
            integrator = new Particle::RK2Integrator<dim, Particle::BaseParticle<dim> >;
          else if (integration_scheme == "rk4")
            integrator = new Particle::RK4Integrator<dim, Particle::BaseParticle<dim> >;
          else if (integration_scheme == "hybrid")
            integrator = new Particle::HybridIntegrator<dim, Particle::BaseParticle<dim> >(&(this->get_triangulation()),
                                                                                           &(this->get_dof_handler()),
                                                                                           &(this->get_mapping()),
                                                                                           &(this->get_solution()));
          else
            Assert (false, ExcNotImplemented());

          // Set up the particle world with the appropriate simulation objects
          world.set_mapping(&(this->get_mapping()));
          world.set_triangulation(&(this->get_triangulation()));
          world.set_dof_handler(&(this->get_dof_handler()));
          world.set_integrator(integrator);
          world.set_mpi_comm(this->get_mpi_communicator());

          // And initialize the world
          world.init();

          next_data_output_time = this->get_time();

          // Add the specified number of particles
          world.global_add_particles(n_initial_tracers);

          initialized = true;
        }

      std::string     result_string = "done", data_file_name;

      // If it's time to generate an output file, call the appropriate functions and reset the timer
      if (this->get_time() >= next_data_output_time)
        {
          set_next_data_output_time (this->get_time());
          data_file_name = output->output_particle_data(world.get_particles(),
                                                        (this->convert_output_to_years() ?
                                                            this->get_time() / year_in_seconds :
                                                            this->get_time()));
          result_string += ". Writing particle graphical output " + data_file_name;
        }

      // Advance the particles in the world by the current timestep
      world.advance_timestep (this->get_timestep(),
                              this->get_solution());

      return std::make_pair("Advecting particles:", result_string);
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
          prm.declare_entry ("Number of tracers", "1000",
                             Patterns::Double (0),
                             "Total number of tracers to create (not per processor or per element). "
                             "The number is parsed as a floating point number (so that one can "
                             "specify, for example, '1e4' particles) but it is interpreted as "
                             "an integer, of course.");
          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "output files. A value of zero indicates that "
                             "output should be generated every time step.\n\n"
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry("Data output format", "none",
                            Patterns::Selection(Particle::Output::output_object_names()),
                            "File format to output raw particle data in.");
          prm.declare_entry("Integration scheme", "rk2",
                            Patterns::Selection("euler|"
                                                "rk2|"
                                                "rk4|"
                                                "hybrid"
                                               ),
                            "Integration scheme to move particles.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    PassiveTracers<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Tracers");
        {
          n_initial_tracers    = static_cast<unsigned int>(prm.get_double ("Number of tracers"));
          data_output_interval = prm.get_double ("Time between data output");
          data_output_format   = prm.get("Data output format");
#ifndef DEAL_II_HAVE_HDF5
          AssertThrow (data_output_format != "hdf5",
                       ExcMessage ("deal.ii was not compiled with HDF5 support, "
                                   "so HDF5 output is not possible. Please "
                                   "recompile deal.ii with HDF5 support turned on."));
#endif
          integration_scheme = prm.get("Integration scheme");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
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
