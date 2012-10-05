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

#include <aspect/global.h>
#include <aspect/postprocess/tracer.h>

namespace aspect
{
    namespace Postprocess
    {
        template <int dim>
        std::pair<std::string,std::string> PassiveTracers<dim>::execute (TableHandler &statistics)
        {
            std::string     result_string = "done.", data_file_name;
            bool            output_data = false;
            
            if (!_initialized)
            {
                // Create an output object depending on what the parameters specify
                if (_data_output_format == "ascii") {
                    _output = new Particle::ASCIIOutput<dim,Particle::BaseParticle<dim> >();
                } else if (_data_output_format == "vtu") {
                    _output = new Particle::VTUOutput<dim,Particle::BaseParticle<dim> >();
                } else if (_data_output_format == "hdf5") {
                    _output = new Particle::HDF5Output<dim,Particle::BaseParticle<dim> >();
                } else {
                    _output = new Particle::NullOutput<dim,Particle::BaseParticle<dim> >();
                }
                
                // Set the output directory for the particle output to be stored in
                _output->set_output_directory(this->get_output_directory());
                
                // Create an integrator object depending on the specified parameter
                if (_integration_scheme == "euler") {
                    _integrator = new Particle::EulerIntegrator<dim, Particle::BaseParticle<dim> >;
                } else if (_integration_scheme == "rk2") {
                    _integrator = new Particle::RK2Integrator<dim, Particle::BaseParticle<dim> >;
                } else if (_integration_scheme == "rk4") {
                    _integrator = new Particle::RK4Integrator<dim, Particle::BaseParticle<dim> >;
                }
                
                // Set up the particle world with the appropriate simulation objects
                _world.set_mapping(&(this->get_mapping()));
                _world.set_triangulation(&(this->get_triangulation()));
                _world.set_dof_handler(&(this->get_dof_handler()));
                _world.set_integrator(_integrator);
                _world.set_mpi_comm(this->get_mpi_communicator());
                
                // And initialize the world
                _world.init();
                
                _next_data_output_time = this->get_time();
                
                // Add the specified number of particles
                _world.global_add_particles(_num_initial_tracers);
                
                _initialized = true;
            }
            
            // Advance the particles in the world by the current timestep
            _world.advance_timestep(this->get_timestep(), this->get_solution());
            
            // If it's time to generate an output file, call the appropriate functions and reset the timer
            if (this->get_time() >= _next_data_output_time)
            {
                set_next_data_output_time (this->get_time());
                data_file_name = _output->output_particle_data(_world.particles(), this->get_time());
                output_data = true;
            }
            if (output_data) result_string += " Wrote particle data: " + data_file_name + ".";
            
            return std::make_pair("Advecting particles...", result_string);
        }
        
        
        
        template <int dim>
        void
        PassiveTracers<dim>::set_next_data_output_time (const double current_time)
        {
            // if output_interval is positive, then set the next output interval to
            // a positive multiple; we need to interpret output_interval either
            // as years or as seconds
            if (_data_output_interval > 0)
            {
                if (this->convert_output_to_years() == true)
                    _next_data_output_time = std::ceil(current_time / (_data_output_interval * year_in_seconds)) *
                    (_data_output_interval * year_in_seconds);
                else
                    _next_data_output_time = std::ceil(current_time / (_data_output_interval)) *
                    (_data_output_interval);
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
                    prm.declare_entry ("Number of tracers", "1e3",
                                       Patterns::Double (0),
                                       "Total number of tracers to create (not per processor or per element).");
                    prm.declare_entry ("Time between data output", "1e8",
                                       Patterns::Double (0),
                                       "The time interval between each generation of "
                                       "output files. A value of zero indicates that "
                                       "output should be generated every time step. "
                                       "Units: years if the "
                                       "'Use years in output instead of seconds' parameter is set; "
                                       "seconds otherwise.");
                    prm.declare_entry("Data output format", "none",
                                      Patterns::Selection("none|"
                                                          "ascii|"
                                                          "vtu|"
                                                          "hdf5"
                                                          ),
                                      "File format to output raw particle data in.");
                    prm.declare_entry("Integration scheme", "rk2",
                                      Patterns::Selection("euler|"
                                                          "rk2|"
                                                          "rk4"
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
                    _num_initial_tracers = prm.get_double ("Number of tracers");
                    _data_output_interval = prm.get_double ("Time between data output");
                    _data_output_format = prm.get("Data output format");
#ifndef DEAL_II_HAVE_HDF5
                    AssertThrow (data_output_format == "hdf5",
                                 ExcMessage ("deal.ii was not compiled with HDF5 support, "
                                             "so HDF5 output is not possible. Please "
                                             "recompile deal.ii with HDF5 support turned on."));
#endif
                    _integration_scheme = prm.get("Integration scheme");
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
