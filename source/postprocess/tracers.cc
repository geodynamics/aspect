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
#include <aspect/postprocess/tracers.h>
#include <aspect/simulator_access.h>

#include <aspect/particle/generator/interface.h>
#include <aspect/particle/integrator/interface.h>
#include <aspect/particle/interpolator/interface.h>
#include <aspect/particle/property/interface.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdio.h>
#include <unistd.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    Tracers<dim>::Tracers ()
      :

      last_output_time(std::numeric_limits<double>::quiet_NaN())
    {}

    template <int dim>
    Tracers<dim>::~Tracers ()
    {}

    template <int dim>
    void
    Tracers<dim>::initialize ()
    {}

    template <int dim>
    void
    Tracers<dim>::generate_and_initialize_particles()
    {
      world.generate_particles();
      world.initialize_particles();
    }

    template <int dim>
    const Particle::World<dim> &
    Tracers<dim>::get_particle_world() const
    {
      return world;
    }

    template <int dim>
    Particle::World<dim> &
    Tracers<dim>::get_particle_world()
    {
      return world;
    }

    template <int dim>
    std::pair<std::string,std::string>
    Tracers<dim>::execute (TableHandler &statistics)
    {
      // if this is the first time we get here, set the last output time
      // to the current time - output_interval. this makes sure we
      // always produce data during the first time step
      if (std::isnan(last_output_time))
        {
          last_output_time = this->get_time() - output_interval;
        }

      // Advance the particles in the world to the current time
      world.advance_timestep();

      statistics.add_value("Number of advected particles",world.n_global_particles());

      // If it's not time to generate an output file or we do not write output
      // return early with the number of particles that were advected
      if ((this->get_time() < last_output_time + output_interval) || !output)
        return std::make_pair("Number of advected particles",Utilities::int_to_string(world.n_global_particles()));


      if (world.get_property_manager().need_update() == Particle::Property::update_output_step)
        world.update_particles();

      set_last_output_time (this->get_time());

      const double output_time = (this->convert_output_to_years() ?
                                  this->get_time() / year_in_seconds :
                                  this->get_time());

      const std::string data_file_name = output->output_particle_data(world.get_particles(),
                                                                      world.get_property_manager().get_data_info(),
                                                                      output_time);

      // record the file base file name in the output file
      statistics.add_value ("Particle file name",
                            this->get_output_directory() + data_file_name);
      return std::make_pair("Writing particle output: ", data_file_name);
    }



    template <int dim>
    void
    Tracers<dim>::set_last_output_time (const double current_time)
    {
      // if output_interval is positive, then update the last supposed output
      // time
      if (output_interval > 0)
        {
          // We need to find the last time output was supposed to be written.
          // this is the last_output_time plus the largest positive multiple
          // of output_intervals that passed since then. We need to handle the
          // edge case where last_output_time+output_interval==current_time,
          // we did an output and std::floor sadly rounds to zero. This is done
          // by forcing std::floor to round 1.0-eps to 1.0.
          const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
          last_output_time = last_output_time + std::floor((current_time-last_output_time)/output_interval*magic) * output_interval/magic;
        }
    }


    template <int dim>
    template <class Archive>
    void Tracers<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &last_output_time
      ;
    }


    template <int dim>
    void
    Tracers<dim>::save (std::map<std::string, std::string> &status_strings) const
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
          if (output)
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
    Tracers<dim>::load (const std::map<std::string, std::string> &status_strings)
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
          if (output)
            output->load(is);
        }
    }


    template <int dim>
    void
    Tracers<dim>::declare_parameters (ParameterHandler &prm)
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
          prm.declare_entry ("Load balancing strategy", "none",
                             Patterns::Selection ("none|remove particles|"
                                                  "remove and add particles|repartition"),
                             "Strategy that is used to balance the computational"
                             "load across processors for adaptive meshes.");
          prm.declare_entry ("Minimum tracers per cell", "100",
                             Patterns::Integer (0),
                             "Limit for how many particles are allowed per cell. This limit is "
                             "useful to prevent coarse cells in adaptive meshes from slowing down "
                             "the whole model. It will only be checked and enforced during "
                             "mesh refinement and MPI transfer of tracers.");
          prm.declare_entry ("Maximum tracers per cell", "100",
                             Patterns::Integer (0),
                             "Limit for how many particles are allowed per cell. This limit is "
                             "useful to prevent coarse cells in adaptive meshes from slowing down "
                             "the whole model. It will only be checked and enforced during "
                             "mesh refinement and MPI transfer of tracers.");
          prm.declare_entry ("Tracer weight", "10",
                             Patterns::Integer (0),
                             "Weight that is associated with the computational load of "
                             "a single particle. The sum of tracer weights will be added "
                             "to the sum of cell weights to determine the partitioning of "
                             "the mesh. Every cell without tracers is associated with a "
                             "weight of 1000.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      Particle::Generator::declare_parameters<dim>(prm);
      Particle::Output::declare_parameters<dim>(prm);
      Particle::Integrator::declare_parameters<dim>(prm);
      Particle::Interpolator::declare_parameters<dim>(prm);
      Particle::Property::Manager<dim>::declare_parameters(prm);
    }


    template <int dim>
    void
    Tracers<dim>::parse_parameters (ParameterHandler &prm)
    {
      // First do some error checking. The current algorithm does not find
      // the cells around particles, if the particles moved more than one
      // cell in one timestep and we are running in parallel, because they
      // skip the layer of ghost cells around our local domain. Assert this
      // is not possible.
      const double CFL_number = prm.get_double ("CFL number");
      const unsigned int n_processes = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

      AssertThrow((n_processes == 1) || (CFL_number <= 1.0),
                  ExcMessage("The current tracer algorithm does not work in "
                             "parallel if the CFL number is larger than 1.0, because "
                             "in this case tracers can move more than one cell's "
                             "diameter in one time step and therefore skip the layer "
                             "of ghost cells around the local subdomain."));

      // Parameters that are handed down to the particle world in this function
      unsigned int max_tracers_per_cell,tracer_weight;
      typename aspect::Particle::World<dim>::ParticleLoadBalancing load_balancing;

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Tracers");
        {
          output_interval = prm.get_double ("Time between data output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;

          max_tracers_per_cell = prm.get_integer("Maximum tracers per cell");
          tracer_weight = prm.get_integer("Tracer weight");

          if (prm.get ("Load balancing strategy") == "none")
            load_balancing = Particle::World<dim>::no_balancing;
          else if (prm.get ("Load balancing strategy") == "remove particles")
            load_balancing = Particle::World<dim>::remove_particles;
          else if (prm.get ("Load balancing strategy") == "remove and add particles")
            load_balancing = Particle::World<dim>::remove_and_add_particles;
          else if (prm.get ("Load balancing strategy") == "repartition")
            load_balancing = Particle::World<dim>::repartition;
          else
            AssertThrow (false, ExcNotImplemented());
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      // Create a generator object depending on what the parameters specify
      Particle::Generator::Interface<dim> *generator
        = Particle::Generator::create_particle_generator<dim> (prm);

      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(generator))
        sim->initialize_simulator (this->get_simulator());
      generator->parse_parameters(prm);


      // Create an output object depending on what the parameters specify
      output.reset(Particle::Output::create_particle_output<dim>
                   (prm));

      // We allow to not generate any output plugin, in which case output is
      // a null pointer. Only initialize output if it was created.
      if (output)
        {
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(output.get()))
            sim->initialize_simulator (this->get_simulator());
          output->parse_parameters(prm);
          output->initialize();
        }


      // Create an integrator object depending on the specified parameter
      Particle::Integrator::Interface<dim> *integrator =
        Particle::Integrator::create_particle_integrator<dim> (prm);

      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(integrator))
        sim->initialize_simulator (this->get_simulator());
      integrator->parse_parameters(prm);

      // Create an interpolator object depending on the specified parameter
      Particle::Interpolator::Interface<dim> *interpolator =
        Particle::Interpolator::create_particle_interpolator<dim> (prm);
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(interpolator))
        sim->initialize_simulator (this->get_simulator());
      interpolator->parse_parameters(prm);


      // Creaty an property_manager object and initialize its properties
      Particle::Property::Manager<dim> *property_manager = new Particle::Property::Manager<dim> ();
      SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(property_manager);
      sim->initialize_simulator (this->get_simulator());
      property_manager->parse_parameters(prm);
      property_manager->initialize();


      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&world))
        sim->initialize_simulator (this->get_simulator());

      // Initialize the particle world with the appropriate settings
      // Ownership of generator, integrator, interpolator and property_manager
      // is transferred to world here
      world.initialize(generator,
                       integrator,
                       interpolator,
                       property_manager,
                       load_balancing,
                       max_tracers_per_cell,
                       tracer_weight);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Tracers,
                                  "tracers",
                                  "A Postprocessor that creates tracer particles that follow the "
                                  "velocity field of the simulation. The particles can be generated "
                                  "and propagated in various ways and they can carry a number of "
                                  "constant or time-varying properties. The postprocessor can write "
                                  "output positions and properties of all tracers at chosen intervals, "
                                  "although this is not mandatory. It also allows other parts of the "
                                  "code to query the tracers for information.")
  }
}
