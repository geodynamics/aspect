/*
  Copyright (C) 2014 - 2024 by the authors of the ASPECT code.

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


#include <aspect/mesh_deformation/external_tool_interface.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    void
    ExternalToolInterface<dim>::
    compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                             AffineConstraints<double> &mesh_velocity_constraints,
                                             const std::set<types::boundary_id> &boundary_ids) const
    {
      // First compute a (global) vector that has the correct velocities
      // set at all boundary nodes:
      const std::vector<std::vector<double>> aspect_surface_velocities = evaluate_aspect_variables_at_points();

      const std::vector<Tensor<1,dim>> external_surface_velocities
        = compute_updated_velocities_at_points(aspect_surface_velocities);

      const LinearAlgebra::Vector v_interpolated
        = interpolate_velocities_to_surface_points(external_surface_velocities);

      // Turn v_interpolated into constraints. For this, loop over all
      // boundary DoFs and if a boundary DoF is locally owned, create a
      // constraint. We later make that consistent across processors to
      // ensure we also know about the locally relevant DoFs'
      // constraints:
      // now insert the relevant part of the solution into the mesh constraints
      const IndexSet constrained_dofs =
        DoFTools::extract_boundary_dofs(mesh_deformation_dof_handler,
                                        ComponentMask(dim, true),
                                        boundary_ids);

      for (const types::global_dof_index index : constrained_dofs)
        {
          if (mesh_velocity_constraints.can_store_line(index))
            if (mesh_velocity_constraints.is_constrained(index)==false)
              {
#if DEAL_II_VERSION_GTE(9,6,0)
                mesh_velocity_constraints.add_constraint(index,
                                                         {},
                                                         v_interpolated(index));
#else
                mesh_velocity_constraints.add_line(index);
                mesh_velocity_constraints.set_inhomogeneity(index, v_interpolated(index));
#endif
              }
        }
    }


    template <int dim>
    void
    ExternalToolInterface<dim>::
    set_evaluation_points (const std::vector<Point<dim>> &evaluation_points)
    {
      // First, save a copy of the points at which we need the solution,
      // among other reasons so that we can track that input arguments
      // for later function calls describe the same number of points.
      this->evaluation_points = evaluation_points;

      // Then also invalidate the previous evaluator:
      remote_point_evaluator.reset();

      // Timo:
      // TODO: Implement set-up phase via MPIRemotePointEvaluation


      // Finally, also ensure that upon mesh refinement, all of the
      // information set herein is invalidated:
      this->get_signals().pre_refinement_store_user_data
      .connect([this](typename parallel::distributed::Triangulation<dim> &)
      {
        this->evaluation_points.clear();
        this->remote_point_evaluator.reset();
      }
              );
    }



    template <int dim>
    std::vector<std::vector<double>>
    ExternalToolInterface<dim>::
    evaluate_aspect_variables_at_points () const
    {
      Assert (remote_point_evaluator != nullptr,
              ExcMessage("You can only call this function if you have previously "
                         "set the evaluation points by calling set_evaluation_points(), "
                         "and if the evaluator has not been invalidated by a mesh "
                         "refinement step."));

      std::vector<std::vector<double>> solution_at_points (evaluation_points.size(),
                                                            std::vector<double>(this->introspection().n_components));
      // Timo: Implement via MPIRemoteEvaluation

      return solution_at_points;
    }



    template <int dim>
    LinearAlgebra::Vector
    ExternalToolInterface<dim>::
    interpolate_velocities_to_surface_points (const std::vector<Tensor<1,dim>> &velocities) const
    {
      Assert (remote_point_evaluator != nullptr,
              ExcMessage("You can only call this function if you have previously "
                         "set the evaluation points by calling set_evaluation_points(), "
                         "and if the evaluator has not been invalidated by a mesh "
                         "refinement step."));
      AssertDimension(velocities.size(), evaluation_points.size());

      // Create the output vector. TODO: We need to get access to the
      // locally owned DoFs index set. Look up how this is done in the
      // other implementations of the Interface base class, if they do
      // it this way at all.

      // LinearAlgebra::Vector vector_with_surface_velocities(this->get_mesh_deformation_handler().mesh_deformation_dof_handler.locally_owned_dofs(),
      //                                             mpi_communicator)
      LinearAlgebra::Vector vector_with_surface_velocities;


      // Timo: Implement via MPIRemoteEvaluation
      (void)velocities;


      return vector_with_surface_velocities;
    }
  }



  namespace MeshDeformation
  {
#define INSTANTIATE(dim) \
  template class ExternalToolInterface<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
