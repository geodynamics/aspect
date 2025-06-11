/*
  Copyright (C) 2018 - 2023 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/openlem.h>

#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    OpenLEM<dim>::OpenLEM()
    {}


    template <int dim>
    void
    OpenLEM<dim>::initialize ()
    {
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          grid_old = Grid<>(10,10);
          grid_new = Grid<>(10,10);

	  // todo: fill the grid new

          // Define all nodes with non-positive elevations (here, the ocean around the
          // island as boundary nodes
          for ( int i = 0; i < grid_new.m; ++i )
            for ( int j = 0; j < grid_new.n; ++j )
              grid_new[i][j].b = grid_new[i][j].h <= 0;

	  // Compute initial flow pattern and water level
	  grid_new.fillLakes();

        }
    }

    template <int dim>
    void
    OpenLEM<dim>::update ()
    {
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years

      // update
      // TODO: get velocity at surface from all processes.
      const double dt = 1e-3; //year?
      for(int i = 0; i < 100; ++i )
      {
        int nc = grid_new.computeFlowDirection();
        printf("Changes in flow direction: %i\n",nc);  
        double ch = grid_new.erode(dt);
        printf("Maximum elevation change: %e\n",ch);
      }
       // Todo: height_diff = grid_new->height - grid_old->height;
       // Todo: interpolate hight_diff to velocties on nodes
       // Todo: move info to all other processes
       // Todo: copy grid_new to grid_old
    }



    /**
     * A function that creates constraints for the velocity of certain mesh
     * vertices (e.g. the surface vertices) for a specific boundary.
     * The calling class will respect
     * these constraints when computing the new vertex positions.
     */
    template <int dim>
    void
    OpenLEM<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                           AffineConstraints<double> &mesh_velocity_constraints,
                                                           const std::set<types::boundary_id> &boundary_ids) const
    {
      // Loop over all boundary indicators to set the velocity constraints
      /*for (const auto boundary_id : boundary_ids)
        VectorTools::interpolate_boundary_values (this->get_mapping(),
                                                  mesh_deformation_dof_handler,
                                                  boundary_id,
                                                  function,
                                                  mesh_velocity_constraints);
      */
    }



    template <int dim>
    bool
    OpenLEM<dim>::
    needs_surface_stabilization () const
    {
      return false;
    }



    template <int dim>
    void OpenLEM<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("openLEM");
        {
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void OpenLEM<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("openLEM");
        {
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(OpenLEM,
                                           "openLEM",
                                           "TODO: A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "Although the function's time variable is interpreted as "
                                           "years when Use years in output instead of seconds is set to true, "
                                           "the boundary deformation velocity should still be given "
                                           "in m/s. The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see {ref}\\`sec:run-aspect:parameters-overview:muparser-format\\`.")
  }
}
