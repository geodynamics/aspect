/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/melt.h>

#include <deal.II/base/mpi.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/fe/mapping_q_cache.h>

#ifdef DEAL_II_WITH_ZLIB
#  include <zlib.h>
#endif

namespace aspect
{
  namespace
  {
    /**
     * Move/rename a file from the given old to the given new name.
     */
    void move_file (const std::string &old_name,
                    const std::string &new_name)
    {
      int error = std::system (("mv " + old_name + " " + new_name).c_str());

      // If the above call failed, e.g. because there is no command-line
      // available, try with internal functions.
      if (error != 0)
        {
          if (Utilities::fexists(new_name))
            {
              error = remove(new_name.c_str());
              AssertThrow (error == 0, ExcMessage(std::string ("Unable to remove file: "
                                                               + new_name
                                                               + ", although it seems to exist. "
                                                               + "The error code is "
                                                               + Utilities::to_string(error) + ".")));
            }

          error = rename(old_name.c_str(),new_name.c_str());
          AssertThrow (error == 0, ExcMessage(std::string ("Unable to rename files: ")
                                              +
                                              old_name + " -> " + new_name
                                              + ". The error code is "
                                              + Utilities::to_string(error) + "."));
        }
    }
  }


  namespace
  {
    /**
     * Save a few of the critical parameters of the current run in the
     * checkpoint file. We will load them again later during
     * restart to verify that they are the same as the ones set
     * in the input file active during restart.
     */
    template <int dim>
    void save_critical_parameters (const Parameters<dim> &parameters,
                                   aspect::oarchive &oa)
    {
      oa << parameters.convert_to_years;
      oa << parameters.surface_pressure;
      oa << parameters.use_operator_splitting;
      oa << parameters.include_melt_transport;
      oa << parameters.stokes_velocity_degree;
      oa << parameters.use_locally_conservative_discretization;
      oa << parameters.use_discontinuous_temperature_discretization;
      oa << parameters.use_discontinuous_composition_discretization;
      oa << parameters.temperature_degree;
      oa << parameters.composition_degrees;
      oa << parameters.pressure_normalization;
      oa << parameters.n_compositional_fields;
      oa << parameters.names_of_compositional_fields;
      oa << parameters.normalized_fields;
      oa << parameters.mesh_deformation_enabled;
      oa << parameters.n_particle_managers;
    }



    /**
     * Load a few of the critical parameters from a checkpoint file during
     * restart to verify that they are the same as the ones currently set
     * in the input file active during restart.
     */
    template <int dim>
    void load_and_check_critical_parameters (const Parameters<dim> &parameters,
                                             aspect::iarchive &ia)
    {
      bool convert_to_years;
      ia >> convert_to_years;
      AssertThrow (convert_to_years == parameters.convert_to_years,
                   ExcMessage ("The value provided for `Use years in output instead of seconds' that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      double surface_pressure;
      ia >> surface_pressure;
      AssertThrow (surface_pressure == parameters.surface_pressure,
                   ExcMessage ("The value of surface pressure that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      bool use_operator_splitting;
      ia >> use_operator_splitting;
      AssertThrow (use_operator_splitting == parameters.use_operator_splitting,
                   ExcMessage ("The operator splitting mode that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      bool include_melt_transport;
      ia >> include_melt_transport;
      AssertThrow (include_melt_transport == parameters.include_melt_transport,
                   ExcMessage ("The melt transport mode that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      unsigned int stokes_velocity_degree;
      ia >> stokes_velocity_degree;
      AssertThrow (stokes_velocity_degree == parameters.stokes_velocity_degree,
                   ExcMessage ("The polynomial degree used for the Stokes "
                               "finite element that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));


      // It is conceivable that one could change this setting from one time
      // step to another, but it is, at best, not tested. So disallow it for
      // now, until someone tests it.
      bool use_locally_conservative_discretization;
      ia >> use_locally_conservative_discretization;
      AssertThrow (use_locally_conservative_discretization == parameters.use_locally_conservative_discretization,
                   ExcMessage ("The value provided for `Use locally conservative discretization' that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      bool use_discontinuous_temperature_discretization;
      ia >> use_discontinuous_temperature_discretization;
      AssertThrow (use_discontinuous_temperature_discretization == parameters.use_discontinuous_temperature_discretization,
                   ExcMessage ("The value provided for `Use discontinuous temperature discretization' that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      std::vector<bool> use_discontinuous_composition_discretization;
      ia >> use_discontinuous_composition_discretization;
      AssertThrow (use_discontinuous_composition_discretization == parameters.use_discontinuous_composition_discretization,
                   ExcMessage ("The value provided for `Use discontinuous composition discretization' that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      unsigned int temperature_degree;
      ia >> temperature_degree;
      AssertThrow (temperature_degree == parameters.temperature_degree,
                   ExcMessage ("The temperature polynomial degree that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      std::vector<unsigned int> composition_degrees;
      ia >> composition_degrees;
      AssertThrow (composition_degrees == parameters.composition_degrees,
                   ExcMessage ("The composition polynomial degree that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      // One could allow changing the pressure normalization between runs, but
      // the change would then lead to a jump in pressure from one time step
      // to the next when we, for example, change from requiring the *surface*
      // average to be zero, to requiring the *domain* average to be zero.
      // That's unlikely what the user really wanted. However, we do allow
      // disabling the pressure normalization in case the user wants to enable
      // a free surface boundary.
      std::string pressure_normalization;
      ia >> pressure_normalization;
      AssertThrow (pressure_normalization == parameters.pressure_normalization ||
                   parameters.pressure_normalization == "no",
                   ExcMessage ("The pressure normalization method that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file and your new "
                               "pressure normalization method is not 'no'. "
                               "The only allowed change for the pressure "
                               "normalization method during a restart is to "
                               "disable normalization."));

      unsigned int n_compositional_fields;
      ia >> n_compositional_fields;
      AssertThrow (n_compositional_fields == parameters.n_compositional_fields,
                   ExcMessage ("The number of compositional fields that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      std::vector<std::string> names_of_compositional_fields;
      ia >> names_of_compositional_fields;
      AssertThrow (names_of_compositional_fields == parameters.names_of_compositional_fields,
                   ExcMessage ("The names of compositional fields that were stored "
                               "in the checkpoint file are not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      std::vector<unsigned int> normalized_fields;
      ia >> normalized_fields;
      AssertThrow (normalized_fields == parameters.normalized_fields,
                   ExcMessage ("The list of normalized fields that was stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      bool mesh_deformation_enabled;
      ia >> mesh_deformation_enabled;
      AssertThrow (mesh_deformation_enabled == parameters.mesh_deformation_enabled,
                   ExcMessage ("The enable mesh deformation settings that were stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));

      unsigned int n_particle_managers;
      ia >> n_particle_managers;
      AssertThrow (n_particle_managers == parameters.n_particle_managers,
                   ExcMessage ("The number of particle systems that were stored "
                               "in the checkpoint file is not the same as the one "
                               "you currently set in your input file. "
                               "These need to be the same during restarting "
                               "from a checkpoint."));
    }
  }


  template <int dim>
  void Simulator<dim>::create_snapshot()
  {
    TimerOutput::Scope timer (computing_timer, "Create snapshot");

    // Take elapsed time from timer so that we can serialize it:
    total_walltime_until_last_snapshot += wall_timer.wall_time();
    wall_timer.restart();

    const unsigned int my_id = Utilities::MPI::this_mpi_process (mpi_communicator);

    // save Triangulation and Solution vectors:
    {
      std::vector<const LinearAlgebra::BlockVector *> x_system
        = { &solution, &old_solution, &old_old_solution };

      // If we are using a deforming mesh, include the mesh velocity, which uses the system dof handler
      if (parameters.mesh_deformation_enabled)
        x_system.push_back( &mesh_deformation->mesh_velocity );

      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
      system_trans (dof_handler);

      system_trans.prepare_for_serialization (x_system);


      // If we are deforming the mesh, also serialize the mesh vertices vector, which
      // uses its own dof handler
      std::vector<const LinearAlgebra::Vector *> x_fs_system;
      std::unique_ptr<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>> mesh_deformation_trans;
      if (parameters.mesh_deformation_enabled)
        {
          x_fs_system.push_back (&mesh_deformation->mesh_displacements);
          x_fs_system.push_back (&mesh_deformation->initial_topography);

          mesh_deformation_trans
            = std::make_unique<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>>
              (mesh_deformation->mesh_deformation_dof_handler);

          mesh_deformation_trans->prepare_for_serialization(x_fs_system);
        }

      signals.pre_checkpoint_store_user_data(triangulation);

      triangulation.save (parameters.output_directory + "restart.mesh.new");
    }

    // save general information This calls the serialization functions on all
    // processes (so that they can take additional action, if necessary, see
    // the manual) but only writes to the restart file on process 0
    {
      std::ostringstream oss;

      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query its properties below.
      {
        aspect::oarchive oa (oss);
        save_critical_parameters (this->parameters, oa);
        oa << (*this);
      }

      // compress with zlib and write to file on the root processor
#ifdef DEAL_II_WITH_ZLIB
      if (my_id == 0)
        {
          uLongf compressed_data_length = compressBound (oss.str().length());
          std::vector<char *> compressed_data (compressed_data_length);
          int err = compress2 (reinterpret_cast<Bytef *>(&compressed_data[0]),
                               &compressed_data_length,
                               reinterpret_cast<const Bytef *>(oss.str().data()),
                               oss.str().length(),
                               Z_BEST_COMPRESSION);
          (void)err;
          Assert (err == Z_OK, ExcInternalError());

          // build compression header
          const uint32_t compression_header[4]
            = { 1,                                   /* number of blocks */
                static_cast<uint32_t>(oss.str().length()), /* size of block */
                static_cast<uint32_t>(oss.str().length()), /* size of last block */
                static_cast<uint32_t>(compressed_data_length)
              }; /* list of compressed sizes of blocks */

          std::ofstream f ((parameters.output_directory + "restart.resume.z.new"));
          f.write(reinterpret_cast<const char *>(compression_header), 4 * sizeof(compression_header[0]));
          f.write(reinterpret_cast<char *>(&compressed_data[0]), compressed_data_length);
          f.close();

          // We check the fail state of the stream _after_ closing the file to
          // make sure the writes were completed correctly. This also catches
          // the cases where the file could not be opened in the first place
          // or one of the write() commands fails, as the fail state is
          // "sticky".
          if (!f)
            AssertThrow(false, ExcMessage ("Writing of the checkpoint file '" + parameters.output_directory
                                           + "restart.resume.z.new' with size "
                                           + Utilities::to_string(4 * sizeof(compression_header[0])+compressed_data_length)
                                           + " failed on processor 0."));
        }
#else
      AssertThrow (false,
                   ExcMessage ("You need to have deal.II configured with the `libz' "
                               "option to support checkpoint/restart, but deal.II "
                               "did not detect its presence when you called `cmake'."));
#endif

    }

    // Wait for everyone to finish writing
    const int ierr = MPI_Barrier(mpi_communicator);
    AssertThrowMPI(ierr);

    // Now rename the snapshots to put the new one in place of the old one.
    // Do this after writing the new one, because writing large checkpoints
    // can be slow, and the model might be cancelled during writing.
    // This way restart remains usable even if restart.new is not completely
    // written.
    if (my_id == 0)
      {
        // if we have previously written a snapshot, then keep the last
        // snapshot in case this one fails to save. Note: static variables
        // will only be initialized once per model run.
        static bool previous_snapshot_exists = (parameters.resume_computation == true);

        if (previous_snapshot_exists == true)
          {
            move_file (parameters.output_directory + "restart.mesh",
                       parameters.output_directory + "restart.mesh.old");
            move_file (parameters.output_directory + "restart.mesh.info",
                       parameters.output_directory + "restart.mesh.info.old");
            move_file (parameters.output_directory + "restart.resume.z",
                       parameters.output_directory + "restart.resume.z.old");

            move_file (parameters.output_directory + "restart.mesh_fixed.data",
                       parameters.output_directory + "restart.mesh_fixed.data.old");

            if (Utilities::fexists(parameters.output_directory + "restart.mesh_variable.data"))
              {
                move_file (parameters.output_directory + "restart.mesh_variable.data",
                           parameters.output_directory + "restart.mesh_variable.data.old");
              }

          }

        move_file (parameters.output_directory + "restart.mesh.new",
                   parameters.output_directory + "restart.mesh");
        move_file (parameters.output_directory + "restart.mesh.new.info",
                   parameters.output_directory + "restart.mesh.info");
        move_file (parameters.output_directory + "restart.resume.z.new",
                   parameters.output_directory + "restart.resume.z");

        move_file (parameters.output_directory + "restart.mesh.new_fixed.data",
                   parameters.output_directory + "restart.mesh_fixed.data");

        if (Utilities::fexists(parameters.output_directory + "restart.mesh.new_variable.data"))
          {
            move_file (parameters.output_directory + "restart.mesh.new_variable.data",
                       parameters.output_directory + "restart.mesh_variable.data");
          }


        // from now on, we know that if we get into this
        // function again that a snapshot has previously
        // been written
        previous_snapshot_exists = true;
      }

    pcout << "*** Snapshot created!" << std::endl << std::endl;
  }



  template <int dim>
  void Simulator<dim>::resume_from_snapshot()
  {
    // By definition, a checkpoint is past the first time step. As a consequence,
    // the Simulator object will not need the initial conditions objects, and
    // we can release the pointers to these objects that we have created in
    // the constructor of this class. If some of the other plugins created there
    // still need access to these initial conditions, they will have created
    // their own shared pointers.
    initial_temperature_manager.reset();
    initial_composition_manager.reset();
#ifdef ASPECT_WITH_WORLD_BUILDER
    // The same applies to the world builder object:
    world_builder.reset();
#endif

    // Then start with the actual deserialization.
    // First check existence of the two restart files
    AssertThrow (Utilities::fexists(parameters.output_directory + "restart.mesh", mpi_communicator),
                 ExcMessage ("You are trying to restart a previous computation, "
                             "but the restart file <"
                             +
                             parameters.output_directory + "restart.mesh"
                             +
                             "> does not appear to exist!"));

    AssertThrow (Utilities::fexists(parameters.output_directory + "restart.resume.z", mpi_communicator),
                 ExcMessage ("You are trying to restart a previous computation, "
                             "but the restart file <"
                             +
                             parameters.output_directory + "restart.resume.z"
                             +
                             "> does not appear to exist!"));

    pcout << "*** Resuming from snapshot!" << std::endl << std::endl;

    // Read resume.z to set up the state of the model
    try
      {
#ifdef DEAL_II_WITH_ZLIB
        const std::string restart_data
          = Utilities::read_and_distribute_file_content (parameters.output_directory + "restart.resume.z",
                                                         mpi_communicator);

        std::istringstream ifs (restart_data);

        uint32_t compression_header[4];
        ifs.read(reinterpret_cast<char *>(compression_header), 4 * sizeof(compression_header[0]));
        Assert(compression_header[0]==1, ExcInternalError());

        std::vector<char> compressed(compression_header[3]);
        std::vector<char> uncompressed(compression_header[1]);
        ifs.read(&compressed[0],compression_header[3]);
        uLongf uncompressed_size = compression_header[1];

        const int err = uncompress(reinterpret_cast<Bytef *>(&uncompressed[0]), &uncompressed_size,
                                   reinterpret_cast<Bytef *>(&compressed[0]), compression_header[3]);
        AssertThrow (err == Z_OK,
                     ExcMessage (std::string("Uncompressing the data buffer resulted in an error with code <")
                                 +
                                 Utilities::int_to_string(err)));

        {
          std::istringstream ss;
          ss.str(std::string (&uncompressed[0], uncompressed_size));

          aspect::iarchive ia (ss);
          load_and_check_critical_parameters(this->parameters, ia);
          ia >> (*this);
        }
#else
        AssertThrow (false,
                     ExcMessage ("You need to have deal.II configured with the `libz' "
                                 "option to support checkpoint/restart, but deal.II "
                                 "did not detect its presence when you called `cmake'."));
#endif
      }
    catch (std::exception &e)
      {
        AssertThrow (false,
                     ExcMessage (std::string("Cannot seem to deserialize the data previously stored!\n")
                                 +
                                 "Some part of the machinery generated an exception that says:\n"
                                 +
                                 e.what()));
      }

    // now that we have resumed from the snapshot load the mesh and solution vectors
    try
      {
        triangulation.load (parameters.output_directory + "restart.mesh");
      }
    catch (...)
      {
        AssertThrow(false, ExcMessage("Cannot open snapshot mesh file or read the triangulation stored there."));
      }

    // if using a cached mapping, update the cache with the new triangulation
    if (MappingQCache<dim> *map = dynamic_cast<MappingQCache<dim>*>(&(*mapping)))
      map->initialize(MappingQGeneric<dim>(4), triangulation);

    setup_dofs();
    global_volume = GridTools::volume (triangulation, *mapping);

    LinearAlgebra::BlockVector distributed_system (system_rhs);
    LinearAlgebra::BlockVector old_distributed_system (system_rhs);
    LinearAlgebra::BlockVector old_old_distributed_system (system_rhs);
    LinearAlgebra::BlockVector distributed_mesh_velocity (system_rhs);

    std::vector<LinearAlgebra::BlockVector *> x_system
      = { &distributed_system, &old_distributed_system, &old_old_distributed_system };

    // If necessary, also include the mesh velocity for deserialization
    // with the system dof handler
    if (parameters.mesh_deformation_enabled)
      x_system.push_back(&distributed_mesh_velocity);

    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
    system_trans (dof_handler);

    system_trans.deserialize (x_system);

    solution = distributed_system;
    old_solution = old_distributed_system;
    old_old_solution = old_old_distributed_system;

    if (parameters.mesh_deformation_enabled)
      {
        // copy the mesh velocity which uses the system dof handler
        mesh_deformation->mesh_velocity = distributed_mesh_velocity;

        // deserialize and copy the vectors using the mesh deformation dof handler
        parallel::distributed::SolutionTransfer<dim, LinearAlgebra::Vector> mesh_deformation_trans( mesh_deformation->mesh_deformation_dof_handler );
        LinearAlgebra::Vector distributed_mesh_displacements( mesh_deformation->mesh_locally_owned,
                                                              mpi_communicator );
        LinearAlgebra::Vector distributed_initial_topography( mesh_deformation->mesh_locally_owned,
                                                              mpi_communicator );
        std::vector<LinearAlgebra::Vector *> fs_system
        = { &distributed_mesh_displacements,
            &distributed_initial_topography
          };

        mesh_deformation_trans.deserialize (fs_system);
        mesh_deformation->mesh_displacements = distributed_mesh_displacements;
        mesh_deformation->initial_topography = distributed_initial_topography;
      }

    signals.post_resume_load_user_data(triangulation);

    // Overwrite the existing statistics file with the one that would have
    // been current at the time of the snapshot we just read back in. We
    // do this because the simulation that created the snapshot may have
    // continued for a few more time steps. The operation here then
    // effectively truncates the 'statistics' file to the position from
    // which the current simulation is going to continue.
    output_statistics();

    // We have to compute the constraints here because the vector that tells
    // us if a cell is a melt cell is not saved between restarts.
    if (parameters.include_melt_transport)
      {
        initialize_current_linearization_point();
        compute_current_constraints ();
      }
  }

}

//why do we need this?!
BOOST_CLASS_TRACKING (aspect::Simulator<2>, boost::serialization::track_never)
BOOST_CLASS_TRACKING (aspect::Simulator<3>, boost::serialization::track_never)


namespace aspect
{

  template <int dim>
  template <class Archive>
  void Simulator<dim>::serialize (Archive &ar, const unsigned int)
  {
    ar &time;
    ar &time_step;
    ar &old_time_step;
    ar &timestep_number;
    ar &pre_refinement_step;
    ar &last_pressure_normalization_adjustment;
    ar &total_walltime_until_last_snapshot;

    ar &postprocess_manager;

    ar &statistics;

    // We do not serialize the statistics_last_write_size and
    // statistics_last_hash variables on purpose. This way, upon
    // restart, they are left at the values initialized by the
    // Simulator::Simulator() constructor, and this causes the
    // Simulator::output_statistics() function to write the
    // whole statistics file anew at the end of the first time
    // step after restart. See there why this is the
    // correct behavior after restart.
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::create_snapshot(); \
  template void Simulator<dim>::resume_from_snapshot();

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
