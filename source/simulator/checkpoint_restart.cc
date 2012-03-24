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


#include <aspect/simulator.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/distributed/solution_transfer.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>



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
    unsigned int my_id = Utilities::System::get_this_mpi_process (mpi_communicator);

    if (my_id == 0)
      {
        // if we have previously written a snapshot, then keep the last
        // snapshot in case this one fails to save
        static bool previous_snapshot_exists = (parameters.resume_computation == true);

        if (previous_snapshot_exists == true)
          {
            move_file (parameters.output_directory + "mesh",
                       parameters.output_directory + "mesh.old");
            move_file (parameters.output_directory + "mesh.info",
                       parameters.output_directory + "mesh.info.old");
            move_file (parameters.output_directory + "resume.txt",
                       parameters.output_directory + "resume.txt.old");

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

      triangulation.save ((parameters.output_directory + "mesh").c_str());
    }

    // save general information
    if (my_id == 0)
      {
        std::ofstream ofs ((parameters.output_directory + "resume.txt").c_str());
        boost::archive::text_oarchive oa (ofs);
        oa << (*this);
      }
    pcout << "*** Snapshot created!" << std::endl;
  }



  template <int dim>
  void Simulator<dim>::resume_from_snapshot()
  {
    triangulation.load ((parameters.output_directory + "mesh").c_str());
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

    std::ifstream ifs ((parameters.output_directory + "resume.txt").c_str());
    boost::archive::text_iarchive ia (ifs);
    ia >> (*this);

    // re-initialize the postprocessors with the current object
    postprocess_manager.initialize (*this);

    pcout << "*** Resuming from snapshot!" << std::endl;
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
  template void Simulator<deal_II_dimension>::create_snapshot();
  template void Simulator<deal_II_dimension>::resume_from_snapshot();
}
