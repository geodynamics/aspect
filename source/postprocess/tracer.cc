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
#include <aspect/simulator_access.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdio.h>
#include <unistd.h>

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
          // Set up the particle world with the appropriate simulation objects
          world.set_integrator(integrator);
          world.set_manager(&property_manager);

          // And initialize the world
          world.init();

          // Let the generator add the specified number of particles if we are
          // not resuming from a snapshot (in that case we have stored particles)
          if (world.get_global_particle_count() == 0)
            {
              generator->generate_particles(world);
              world.initialize_particles();
            }

          initialized = true;
          next_data_output_time = this->get_time();
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
    template <class Archive>
    void PassiveTracers<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &next_data_output_time
      ;

      // We do not serialize mesh_changed but use the default (true) from our
      // constructor. This will result in a new mesh file the first time we
      // create visualization output after resuming from a snapshot. Otherwise
      // we might get corrupted graphical output, because the ordering of
      // vertices can be different after resuming.
    }


    template <int dim>
    void
    PassiveTracers<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;
      aspect::oarchive oa (os);

      oa << world;

      const unsigned int mpi_tag = 124;

      // on processor 0, collect all of the data the individual processors send
       // and concatenate them as serialized strings into the archive
       if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
         {
           // loop through all of the other processors and collect
           // data, then write it in the archive
           // TODO: this assumes that the number of processes does not change
           for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
             {
               // get the data length and data
               MPI_Status status;
               MPI_Probe(p, mpi_tag, this->get_mpi_communicator(), &status);

               int data_length;
               MPI_Get_count(&status, MPI_CHAR, &data_length);

               std::string tmp(data_length,'\0');
               MPI_Recv (&tmp[0], data_length, MPI_CHAR, p, mpi_tag,
                         this->get_mpi_communicator(), &status);

               oa << tmp;
             }

           oa << (*this);
           output->save(os);
         }
       else
         // on other processors, send the serialized data to processor zero.
         {
           MPI_Send (&os.str()[0], os.str().size(), MPI_CHAR, 0, mpi_tag,
                     this->get_mpi_communicator());
         }

      status_strings["Tracers"] = os.str();

    }


    template <int dim>
    void
    PassiveTracers<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("Tracers") != status_strings.end())
        {
          std::istringstream is (status_strings.find("Tracers")->second);
          aspect::iarchive ia (is);

          // Load the particle world of the first process. This will be
          // overwritten on all processes but the first one later on, but every
          // process needs to load the whole archive to correctly deserialize
          // the data
          ia >> world;

          // Then loop through all of the other processors, but save only
          // the data corresponding to this process. The tracers might not
          // be in the domain of this process anymore due to a freshly build
          // mesh, but this does not matter, because they will simply be send
          // around until every tracer has found its process in the usual way
          // after a mesh change.
          // TODO: this assumes that the number of processes does not change
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              std::string tmp;
              ia >> tmp;

              if (p == Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                {
                  std::istringstream ws (tmp);
                  aspect::iarchive wa (ws);
                  wa >> world;
                }
            }

          ia >> (*this);
          output->load(is);
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
      output->initialize();

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
