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


#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <deal.II/base/mpi.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/solution_transfer.h>

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
      const int error = system (("mv " + old_name + " " + new_name).c_str());

      AssertThrow (error == 0, ExcMessage(std::string ("Can't move files: ")
                                          +
                                          old_name + " -> " + new_name));
    }
  }

    template<int dim>
    void Simulator<dim>::save_triangulation(const std::string file_name)
    {
      // save Triangulation and Solution vectors:
      {
        std::vector<const LinearAlgebra::BlockVector *> x_system (3);
        x_system[0] = &solution;
        x_system[1] = &old_solution;
        x_system[2] = &old_old_solution;

        //If we are using a free surface, include the mesh velocity, which uses the system dof handler
        if (parameters.free_surface_enabled)
          x_system.push_back( &free_surface->mesh_velocity );

        parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
                system_trans (dof_handler);

        system_trans.prepare_serialization (x_system);

        //If we are using a free surface, also serialize the mesh vertices vector, which
        //uses its own dof handler
        std::vector<const LinearAlgebra::Vector *> x_fs_system (2);
        std_cxx11::unique_ptr<parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector> > freesurface_trans;
        if (parameters.free_surface_enabled)
        {
          freesurface_trans.reset (new parallel::distributed::SolutionTransfer<dim,LinearAlgebra::Vector>
                                           (free_surface->free_surface_dof_handler));

          x_fs_system[0] = &free_surface->mesh_vertices;
          x_fs_system[1] = &free_surface->mesh_vertex_velocity;

          freesurface_trans->prepare_serialization(x_fs_system);
        }

        triangulation.save ((parameters.output_directory + file_name).c_str());
      }
    }

    template<int dim>
    void Simulator<dim>::serialize_all(const std::string file_name)
    {
      unsigned int my_id = dealii::Utilities::MPI::this_mpi_process (mpi_communicator);
      // save general information This calls the serialization functions on all
      // processes (so that they can take additional action, if necessary, see
      // the manual) but only writes to the restart file on process 0
      {
        std::ostringstream oss;

        // serialize into a stringstream
        aspect::oarchive oa (oss);
        oa << (*this);

        // compress with zlib and write to file on the root processor
#ifdef DEAL_II_WITH_ZLIB
        if (my_id == 0)
        {
          uLongf compressed_data_length = compressBound (oss.str().length());
          std::vector<char *> compressed_data (compressed_data_length);
          int err = compress2 ((Bytef *) &compressed_data[0],
                               &compressed_data_length,
                               (const Bytef *) oss.str().data(),
                               oss.str().length(),
                               Z_BEST_COMPRESSION);
          (void)err;
          Assert (err == Z_OK, ExcInternalError());

          // build compression header
          const uint32_t compression_header[4]
                  = { 1,                                   /* number of blocks */
                      (uint32_t)oss.str().length(), /* size of block */
                      (uint32_t)oss.str().length(), /* size of last block */
                      (uint32_t)compressed_data_length
                  }; /* list of compressed sizes of blocks */

          std::ofstream f ((parameters.output_directory + file_name).c_str());
          f.write((const char *)compression_header, 4 * sizeof(compression_header[0]));
          f.write((char *)&compressed_data[0], compressed_data_length);
        }
#else
        AssertThrow (false,
                   ExcMessage ("You need to have deal.II configured with the 'libz' "
                               "option to support checkpoint/restart, but deal.II "
                               "did not detect its presence when you called 'cmake'."));
#endif

      }
    }

    template<int dim>
    void Simulator<dim>::log_checkpoint(const std::string filename_for_triangulation, const std::string filename_for_serialization, bool is_quicksave)
    {
      std::ofstream checkpoint_log;
      std::string checkpoint_file_name = parameters.output_directory + checkpoint_log_file;

      if(!Utilities::fexists(checkpoint_file_name))
      {
        checkpoint_log.open(checkpoint_file_name, std::ios_base::out);
        checkpoint_log << "Time_step_number filename_triangulation filename_mesh quicksave_slot_id" << std::endl;
      }
      else
        checkpoint_log.open(checkpoint_file_name, std::ios_base::app);

      if (is_quicksave)
        checkpoint_log << timestep_number << " "
                       << filename_for_triangulation << " "
                       << filename_for_serialization << " "
                       << dealii::Utilities::int_to_string(n_quicksaves % parameters.quicksave_slots) << std::endl;
      else
        checkpoint_log << timestep_number << " "
        << filename_for_triangulation << " "
        << filename_for_serialization << " "
        << "-" << std::endl;

      checkpoint_log.close();
    }

    template <int dim>
    void Simulator<dim>::create_snapshot()
    {
      computing_timer.enter_section ("Create snapshot");

      std::string filename_for_triangulation = "restart.mesh-" + dealii::Utilities::int_to_string(timestep_number);
      std::string filename_for_serialization = "restart.resume-" + dealii::Utilities::int_to_string(timestep_number) + ".z";

      save_triangulation(filename_for_triangulation);
      serialize_all(filename_for_serialization);

      log_checkpoint(filename_for_triangulation, filename_for_serialization, false);

      pcout << "*** Snapshot created!" << std::endl << std::endl;
      computing_timer.exit_section();
    }

    template <int dim>
    void Simulator<dim>::quicksave_snapshot()
    {
      computing_timer.enter_section ("Create snapshot");

      std::string filename_for_triangulation = "quicksave-slot-" + dealii::Utilities::int_to_string(n_quicksaves % parameters.quicksave_slots);
      std::string filename_for_serialization = "quicksave-slot-" + dealii::Utilities::int_to_string(n_quicksaves % parameters.quicksave_slots) + ".z";

      save_triangulation(filename_for_triangulation);
      serialize_all(filename_for_serialization);

      log_checkpoint(filename_for_triangulation, filename_for_serialization, true);

      n_quicksaves++;

      pcout << "*** Quickly saving snapshot! Snapshot Created!" << std::endl << std::endl;
      computing_timer.exit_section();
    }


  template <int dim>
  void Simulator<dim>::resume_from_snapshot()
  {
    std::string filename_for_triangulation;
    std::string filename_to_deserialize;

    if (parameters.resume_from_tsn != 0) {
      //Check for checkpoint files in output directory
      if (Utilities::fexists(parameters.output_directory + "restart.mesh-" + dealii::Utilities::to_string(parameters.resume_from_tsn)) &&
              Utilities::fexists(parameters.output_directory + "restart.resume-" + dealii::Utilities::to_string(parameters.resume_from_tsn) + ".z")) {
        filename_for_triangulation = parameters.output_directory + "restart.mesh-" + dealii::Utilities::to_string(parameters.resume_from_tsn);
        filename_to_deserialize = parameters.output_directory + "restart.resume-" + dealii::Utilities::to_string(parameters.resume_from_tsn) + ".z";
      }
        //Check for checkpoint files in run directory
      else if (Utilities::fexists("restart.mesh-" + dealii::Utilities::to_string(parameters.resume_from_tsn)) &&
              Utilities::fexists("restart.resume-" + dealii::Utilities::to_string(parameters.resume_from_tsn) + ".z")) {
        filename_for_triangulation = "restart.mesh-" + dealii::Utilities::to_string(parameters.resume_from_tsn);
        filename_to_deserialize = "restart.resume-" + dealii::Utilities::to_string(parameters.resume_from_tsn) + ".z";
      }
    }
      // Try to pick up the latest checkpoint file from checkpoint.log in the output directory.
    else if (Utilities::fexists(parameters.output_directory + "checkpoint.log"))
    {
      std::ifstream checkpoint_log(parameters.output_directory + "checkpoint.log", std::ios_base::in);
      if (!checkpoint_log) {
        AssertThrow (false,
                     ExcMessage("Failed to find the checkpoint log file in " + parameters.output_directory +
                                ". Please rerun aspect or specify a valid file name to resume from."
                     ))
      }

      std::string last_line = Utilities::get_last_line(&checkpoint_log);
      std::vector<std::string> tokens;
      std::string delimiter = " ";

      size_t pos = 0;
      while ((pos = last_line.find(delimiter)) != std::string::npos) {
        tokens.push_back(last_line.substr(0, pos));
        last_line.erase(0, pos + delimiter.length());
      }

      filename_for_triangulation = parameters.output_directory + tokens[1];
      filename_to_deserialize = parameters.output_directory + tokens[2];
    }
      // For backward compatability purposes, we restart from the below checkpoint files.
    else
    {
      filename_for_triangulation = parameters.output_directory + "restart.mesh";
      filename_to_deserialize = parameters.output_directory + "restart.resume.z";
    }


    {
      std::ifstream in(filename_for_triangulation.c_str());
      if (!in) AssertThrow (false,
                            ExcMessage(std::string("You are trying to restart a previous computation, "
                                                           "but the restart file <")
                                       +
                                       filename_for_triangulation
                                       +
                                       "> does not appear to exist!"));
    }

    {
      std::ifstream in (filename_to_deserialize.c_str());
      if (!in)
      AssertThrow (false,
                   ExcMessage (std::string("You are trying to restart a previous computation, "
                                                   "but the restart file <")
                               +
                               filename_to_deserialize
                               +
                               "> does not appear to exist!"));
    }

    pcout << "*** Resuming from snapshot!" << std::endl << std::endl;

    try
    {
      triangulation.load ((filename_for_triangulation).c_str());
    }
    catch (...)
    {
      AssertThrow(false, ExcMessage("Cannot open snapshot mesh file or read the triangulation stored there."));
    }
    global_volume = GridTools::volume (triangulation, mapping);
    setup_dofs();

    LinearAlgebra::BlockVector
            distributed_system (system_rhs);
    LinearAlgebra::BlockVector
            old_distributed_system (system_rhs);
    LinearAlgebra::BlockVector
            old_old_distributed_system (system_rhs);
    LinearAlgebra::BlockVector
            distributed_mesh_velocity (system_rhs);

    std::vector<LinearAlgebra::BlockVector *> x_system (3);
    x_system[0] = & (distributed_system);
    x_system[1] = & (old_distributed_system);
    x_system[2] = & (old_old_distributed_system);
    //If necessary, also include the mesh velocity for deserialization
    //with the system dof handler
    if (parameters.free_surface_enabled)
      x_system.push_back(&distributed_mesh_velocity);

    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
            system_trans (dof_handler);

    system_trans.deserialize (x_system);

    solution = distributed_system;
    old_solution = old_distributed_system;
    old_old_solution = old_old_distributed_system;

    if (parameters.free_surface_enabled)
    {
      //copy the mesh velocity which uses the system dof handler
      free_surface->mesh_velocity = distributed_mesh_velocity;

      //deserialize and copy the vectors using the free surface dof handler
      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::Vector> freesurface_trans( free_surface->free_surface_dof_handler );
      LinearAlgebra::Vector distributed_mesh_vertices( free_surface->mesh_locally_owned,
                                                       mpi_communicator );
      LinearAlgebra::Vector distributed_mesh_vertex_velocity( free_surface->mesh_locally_owned,
                                                              mpi_communicator );
      std::vector<LinearAlgebra::Vector *> fs_system(2);
      fs_system[0] = &distributed_mesh_vertices;
      fs_system[1] = &distributed_mesh_vertex_velocity;

      freesurface_trans.deserialize (fs_system);
      free_surface->mesh_vertices = distributed_mesh_vertices;
      free_surface->mesh_vertex_velocity = distributed_mesh_vertex_velocity;
      //Make sure the mesh conforms to the mesh_vertices vector
      free_surface->displace_mesh();
      free_surface->detach_manifolds();
    }


    // read zlib compressed resume.z
    try
    {
#ifdef DEAL_II_WITH_ZLIB
      std::ifstream ifs ((filename_to_deserialize).c_str());
      AssertThrow(ifs.is_open(),
                  ExcMessage("Cannot open snapshot resume file."));

      uint32_t compression_header[4];
      ifs.read((char *)compression_header, 4 * sizeof(compression_header[0]));
      Assert(compression_header[0]==1, ExcInternalError());

      std::vector<char> compressed(compression_header[3]);
      std::vector<char> uncompressed(compression_header[1]);
      ifs.read(&compressed[0],compression_header[3]);
      uLongf uncompressed_size = compression_header[1];

      const int err = uncompress((Bytef *)&uncompressed[0], &uncompressed_size,
                                 (Bytef *)&compressed[0], compression_header[3]);
      AssertThrow (err == Z_OK,
                   ExcMessage (std::string("Uncompressing the data buffer resulted in an error with code <")
                               +
                               dealii::Utilities::int_to_string(err)));

      {
        std::istringstream ss;
        ss.str(std::string (&uncompressed[0], uncompressed_size));
        aspect::iarchive ia (ss);

        ia >> (*this);
      }
#else
      AssertThrow (false,
                     ExcMessage ("You need to have deal.II configured with the 'libz' "
                                 "option to support checkpoint/restart, but deal.II "
                                 "did not detect its presence when you called 'cmake'."));
#endif
    }
    catch (std::exception &e)
    {
      AssertThrow (false,
                   ExcMessage (std::string("Cannot seem to deserialize the data previously stored!\n")
                               +
                               "Some part of the machinery generated an exception that says <"
                               +
                               e.what()
                               +
                               ">"));
    }
  }
}

//why do we need this?!
BOOST_CLASS_TRACKING (aspect::Simulator<2>, boost::serialization::track_never)
BOOST_CLASS_TRACKING (aspect::Simulator<3>, boost::serialization::track_never)


namespace aspect
{

  template <int dim>
  template<class Archive>
  void Simulator<dim>::serialize (Archive &ar, const unsigned int)
  {
    ar &time;
    ar &time_step;
    ar &old_time_step;
    ar &timestep_number;
    ar &pre_refinement_step;
    ar &pressure_adjustment;
    ar &n_quicksaves;

    ar &postprocess_manager &statistics;

  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::create_snapshot(); \
  template void Simulator<dim>::quicksave_snapshot(); \
  template void Simulator<dim>::resume_from_snapshot();

  ASPECT_INSTANTIATE(INSTANTIATE)
}
