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

#include <deal.II/base/mpi.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/solution_transfer.h>

#include <zlib.h>

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


  template <int dim>
  void Simulator<dim>::create_snapshot()
  {
    computing_timer.enter_section ("Create snapshot");
    unsigned int my_id = Utilities::MPI::this_mpi_process (mpi_communicator);

    if (my_id == 0)
      {
        // if we have previously written a snapshot, then keep the last
        // snapshot in case this one fails to save
        static bool previous_snapshot_exists = (parameters.resume_computation == true);

        if (previous_snapshot_exists == true)
          {
            move_file (parameters.output_directory + "restart.mesh",
                       parameters.output_directory + "restart.mesh.old");
            move_file (parameters.output_directory + "restart.mesh.info",
                       parameters.output_directory + "restart.mesh.info.old");
            move_file (parameters.output_directory + "restart.resume.z",
                       parameters.output_directory + "restart.resume.z.old");

            // from now on, we know that if we get into this
            // function again that a snapshot has previously
            // been written
            previous_snapshot_exists = true;
          }
      }

    // save Triangulation and Solution vectors:
    {
      std::vector<const LinearAlgebra::BlockVector *> x_system (3);
      x_system[0] = &solution;
      x_system[1] = &old_solution;
      x_system[2] = &old_old_solution;

      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
      system_trans (dof_handler);

      system_trans.prepare_serialization (x_system);

      triangulation.save ((parameters.output_directory + "restart.mesh").c_str());
    }

    // save general information This calls the serialization functions on all
    // processes (so that they can take additional action, if necessary, see
    // the manual) but only writes to the restart file on process 0
    {
      std::ostringstream oss;

      // serialize into a stringstream
      aspect::oarchive oa (oss);
      oa << (*this);

      // compress with zlib and write to file on the root processor
      if (my_id == 0)
        {
          uLongf compressed_data_length = compressBound (oss.str().length());
          std::vector<char *> compressed_data (compressed_data_length);
          int err = compress2 ((Bytef *) &compressed_data[0],
                               &compressed_data_length,
                               (const Bytef *) oss.str().data(),
                               oss.str().length(),
                               Z_BEST_COMPRESSION);
          Assert (err == Z_OK, ExcInternalError());

          // build compression header
          const uint32_t compression_header[4]
            = { 1,                                   /* number of blocks */
                (uint32_t)oss.str().length(), /* size of block */
                (uint32_t)oss.str().length(), /* size of last block */
                (uint32_t)compressed_data_length
              }; /* list of compressed sizes of blocks */

          std::ofstream f ((parameters.output_directory + "restart.resume.z").c_str());
          f.write((const char *)compression_header, 4 * sizeof(compression_header[0]));
          f.write((char *)&compressed_data[0], compressed_data_length);
        }

    }
    pcout << "*** Snapshot created!" << std::endl << std::endl;
    computing_timer.exit_section();
  }



  template <int dim>
  void Simulator<dim>::resume_from_snapshot()
  {
    // first check existence of the two restart files
    {
      const std::string filename = parameters.output_directory + "restart.mesh";
      std::ifstream in (filename.c_str());
      if (!in)
        AssertThrow (false,
                     ExcMessage (std::string("You are trying to restart a previous computation, "
                                             "but the restart file <")
                                 +
                                 filename
                                 +
                                 "> does not appear to exist!"));
    }
    {
      const std::string filename = parameters.output_directory + "restart.resume.z";
      std::ifstream in (filename.c_str());
      if (!in)
        AssertThrow (false,
                     ExcMessage (std::string("You are trying to restart a previous computation, "
                                             "but the restart file <")
                                 +
                                 filename
                                 +
                                 "> does not appear to exist!"));
    }

    pcout << "*** Resuming from snapshot!" << std::endl << std::endl;

    try
      {
        triangulation.load ((parameters.output_directory + "restart.mesh").c_str());
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
    std::vector<LinearAlgebra::BlockVector *> x_system (3);
    x_system[0] = & (distributed_system);
    x_system[1] = & (old_distributed_system);
    x_system[2] = & (old_old_distributed_system);

    parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector>
    system_trans (dof_handler);

    system_trans.deserialize (x_system);

    solution = distributed_system;
    old_solution = old_distributed_system;
    old_old_solution = old_old_distributed_system;


    // read zlib compressed resume.z
    try
      {
        std::ifstream ifs ((parameters.output_directory + "restart.resume.z").c_str());
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
                                 Utilities::int_to_string(err)));

        {
          std::istringstream ss;
          ss.str(std::string (&uncompressed[0], uncompressed_size));
          aspect::iarchive ia (ss);
          ia >> (*this);
        }
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

    ar &postprocess_manager &statistics;

  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::create_snapshot(); \
  template void Simulator<dim>::resume_from_snapshot();

  ASPECT_INSTANTIATE(INSTANTIATE)
}
