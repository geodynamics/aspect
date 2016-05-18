/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
#include <aspect/prescribed_stokes_solution/binary_data.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/dofs/dof_renumbering.h>

namespace aspect
{
  namespace PrescribedStokesSolution
  {
    template <int dim>
    BinaryData<dim>::BinaryData ()
    {}


    template <int dim>
    void
    BinaryData<dim>::initialize ()
    {
      triangulation = new parallel::distributed::Triangulation<dim>(this->get_mpi_communicator()); 
    }

    template <int dim>
    void
    BinaryData<dim>::update ()
    {
      std::string fileName = dataDirectory + "/" + "solution-" + Utilities::int_to_string(this->get_timestep_number(), 5) + ".mesh";
      triangulation->copy_triangulation(this->get_triangulation());
      triangulation->load(fileName.c_str(), true);
      DoFHandler<dim> dof_handler (*triangulation);
      dof_handler.distribute_dofs(this->get_fe());
      DoFRenumbering::hierarchical (dof_handler);
      DoFRenumbering::component_wise (dof_handler,
                                      this->introspection().get_components_to_blocks());
      LinearAlgebra::BlockVector tmp_solution (this->introspection().index_sets.system_partitioning, this->get_mpi_communicator());
      parallel::distributed::SolutionTransfer<dim, LinearAlgebra::BlockVector> sol_trans(dof_handler);
      sol_trans.deserialize(tmp_solution);

      std::string fileNameTmp = this->get_output_directory() + "/" + "solution-" + Utilities::int_to_string(this->get_timestep_number(), 5) + ".txt"; 
      std::ofstream output(fileNameTmp);
      dealii::BlockVector<double> solTmp(tmp_solution);
      solTmp.print(output);
      output.close(); 

      triangulation->clear();
    }

    template <int dim>
    void
    BinaryData<dim>::stokes_solution (const Point<dim> &position, Vector<double> &value) const
    {
      value(0) = 0; 
      value(1) = 0; 
      if (dim == 3)
        value(2) = 0; 
      value(dim) = 0;  // makes pressure 0, must set pressure
    }


    template <int dim>
    void
    BinaryData<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
      {
        prm.enter_subsection("Binary data model");
        {
          prm.declare_entry("Data directory", "$ASPECT_SOURCE_DIR/data/prescribed-stokes-solution/", Patterns::DirectoryName (),
                                                          "Location of binary data on local filesystem.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    BinaryData<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Prescribed Stokes solution");
      {
        prm.enter_subsection("Binary data model");
        {
          dataDirectory = prm.get("Data directory");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedStokesSolution
  {
    ASPECT_REGISTER_PRESCRIBED_STOKES_SOLUTION(BinaryData,
                                               "binary data",
                                               "Implementation of a model in which the velocity"
                                               "is derived from files containing data "
                                               "in binary format. Note the required format of the "
                                               "input data: The first lines may contain any number of comments "
                                               "if they begin with '#', but one of these lines needs to "
                                               "contain the number of grid points in each dimension as "
                                               "for example '# POINTS: 3 3'. "
                                               "The order of the data columns "
                                               "has to be 'x', 'y', 'v$_x$' , 'v$_y$' in a 2d model and "
                                               " 'x', 'y', 'z', 'v$_x$' , 'v$_y$' , 'v$_z$' in a 3d model. "
                                               "Note that the data in the input "
                                               "files need to be sorted in a specific order: "
                                               "the first coordinate needs to ascend first, "
                                               "followed by the second and the third at last in order to "
                                               "assign the correct data to the prescribed coordinates. "
                                               "If you use a spherical model, "
                                               "then the data will still be handled as Cartesian, "
                                               "however the assumed grid changes. 'x' will be replaced by "
                                               "the radial distance of the point to the bottom of the model, "
                                               "'y' by the azimuth angle and 'z' by the polar angle measured "
                                               "positive from the north pole. The grid will be assumed to be "
                                               "a latitude-longitude grid. Note that the order "
                                               "of spherical coordinates is 'r', 'phi', 'theta' "
                                               "and not 'r', 'theta', 'phi', since this allows "
                                               "for dimension independent expressions.")
  }
}
