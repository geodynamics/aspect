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
#include <aspect/simulator_signals.h>

#include <deal.II/numerics/vector_tools_evaluate.h>

/*
TODO:
- initial topography from external tool: for now compute_initial_deformation_on_boundary(), later refactor
  to provide FE vector
*/

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
        = interpolate_external_velocities_to_surface_support_points(external_surface_velocities);

      // TODO: need ghost values of v_interpolated?

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

      // TODO: make consistent?
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

      // Set up RemotePointEvaluation. The evaluation points are given in reference coordinates,
      // so we need to use a simple mapping instead of the one stored in the Simulator. The latter
      // would produce the deformed mesh. We currently always use a Q1 mapping when mesh deformation
      // is enabled, so a Q1 mapping is the right choice.
      static MappingQ<dim> mapping(1);
      remote_point_evaluator = std::make_unique<Utilities::MPI::RemotePointEvaluation<dim, dim>>();
      remote_point_evaluator->reinit(this->evaluation_points, this->get_triangulation(), mapping);


      // Create a mapping from evaluation points to support points. Note that one evaluation point can map to
      // multiple support points.
      {
        // For now, we only support a single MPI rank for ASPECT and the external mesh deformation class. Because deciding
        // which evaluation point is closest to a support point requires somewhat complex MPI communication.
        Assert(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
               ExcNotImplemented("This function is not implemented for parallel computations."));

        const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();

        std::vector<double> squared_distances(mesh_dof_handler.locally_owned_dofs().size(), std::numeric_limits<double>::max());
        std::vector<DofToEvalPointData> closest_evaluation_point_and_component(mesh_dof_handler.locally_owned_dofs().size(),
                                                                               DofToEvalPointData {numbers::invalid_dof_index, numbers::invalid_unsigned_int, numbers::invalid_unsigned_int});

        // TODO: do we need to support the case of more than one different mesh deformation plugin to be active?
        const auto boundary_ids = this->get_mesh_deformation_boundary_indicators();

        const IndexSet boundary_dofs = DoFTools::extract_boundary_dofs(mesh_dof_handler, ComponentMask(dim, true), boundary_ids);

        const unsigned int dofs_per_cell = mesh_dof_handler.get_fe().dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

        // The remote_point_evaluator will gives us the velocities in all evaluation points that are within one of our locally
        // owned cells. The lambda defined below receives a list of points and their velocities for each cell. The coordinates
        // are given in coordinates of the unit cell.
        // For each support point of the velocity DoFHandler, we will try to find the closest evaluation point. We
        // do this by keeping track of the squared distance of the closest evaluation point checked so far.


        // Note: We assume that process_and_evaluate() does not call our lambda concurrently, otherwise we would have write
        // conflicts when updating closest_evaluation_point_and_component and squared_distances.

        const auto eval_func = [&](const ArrayView<const unsigned int> &values,
                                   const typename Utilities::MPI::RemotePointEvaluation<dim>::CellData &cell_data)
        {
          std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
          for (const auto cell_index : cell_data.cell_indices())
            {
              const auto cell_dofs =
                cell_data.get_active_cell_iterator(cell_index)->as_dof_handler_iterator(
                  mesh_dof_handler);
              cell_dofs->get_dof_indices(cell_dof_indices);

              const ArrayView<const Point<dim>> unit_points = cell_data.get_unit_points(cell_index);
              const auto local_values = cell_data.get_data_view(cell_index, values);

              const std::vector< Point< dim >> &support_points = mesh_dof_handler.get_fe().get_unit_support_points();
              for (unsigned int i=0; i<unit_points.size(); ++i)
                {
                  for (unsigned int j=0; j<support_points.size(); ++j)
                    {
                      const double distance_sq = unit_points[i].distance_square(support_points[j]);
                      if (distance_sq < squared_distances[cell_dof_indices[j]])
                        {
                          squared_distances[cell_dof_indices[j]] = distance_sq;
                          const unsigned int component = mesh_dof_handler.get_fe().system_to_component_index(j).first;
                          closest_evaluation_point_and_component[cell_dof_indices[j]] =
                            DofToEvalPointData {cell_dof_indices[j], local_values[i], component};
                        }
                    }
                }
            }
        };

        std::vector<unsigned int> indices (evaluation_points.size());
        std::iota(indices.begin(), indices.end(), 0);
        this->remote_point_evaluator->template process_and_evaluate<unsigned int, 1>(indices, eval_func, /*sort_data*/ true);

        map_dof_to_eval_point.clear();
        for (unsigned int i=0; i<closest_evaluation_point_and_component.size(); ++i)
          {
            if (closest_evaluation_point_and_component[i].dof_index != numbers::invalid_dof_index)
              map_dof_to_eval_point.push_back(closest_evaluation_point_and_component[i]);
          }
      }

      // Finally, also ensure that upon mesh refinement, all of the
      // information set herein is invalidated:
      this->get_signals().pre_refinement_store_user_data
      .connect([this](typename parallel::distributed::Triangulation<dim> &)
      {
        this->evaluation_points.clear();
        this->remote_point_evaluator.reset();
      });
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

      const unsigned int n_components = this->introspection().n_components;
      std::vector<std::vector<double>> solution_at_points (evaluation_points.size(), std::vector<double>(n_components, 0.0));

      // VectorTools::point_values can evaluate N components at a time, but this is a template argument and not a
      // runtime argument. For now, we just evaluate them one component at the time. Of course it would be more
      // efficient to branch and evaluate up to K at a time (for a reasonable number of K, say 10). Maybe something
      // to put directly into deal.II...
      for (unsigned int c=0; c<n_components; ++c)
        {
          const std::vector<double> values = VectorTools::point_values<1>(*this->remote_point_evaluator,
                                                                          this->get_dof_handler(),
                                                                          this->get_solution(),
                                                                          dealii::VectorTools::EvaluationFlags::avg,
                                                                          c);
          for (unsigned int i=0; i<evaluation_points.size(); ++i)
            solution_at_points[i][c] = values[i];
        }

      return solution_at_points;
    }



    template <int dim>
    LinearAlgebra::Vector
    ExternalToolInterface<dim>::
    interpolate_external_velocities_to_surface_support_points (const std::vector<Tensor<1,dim>> &velocities) const
    {
      Assert (remote_point_evaluator != nullptr,
              ExcMessage("You can only call this function if you have previously "
                         "set the evaluation points by calling set_evaluation_points(), "
                         "and if the evaluator has not been invalidated by a mesh "
                         "refinement step."));
      AssertDimension(velocities.size(), evaluation_points.size());


      // Create the output vector.
      const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();
      LinearAlgebra::Vector vector_with_surface_velocities(mesh_dof_handler.locally_owned_dofs(),
                                                           this->get_mpi_communicator());

      for (const auto &entry : map_dof_to_eval_point)
        vector_with_surface_velocities[entry.dof_index] = velocities[entry.evaluation_point_index][entry.component];

      vector_with_surface_velocities.compress(VectorOperation::insert);

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
