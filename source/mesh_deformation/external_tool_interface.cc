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
- PointDataOut: move to deal.II
- initial topography from external tool: for now compute_initial_deformation_on_boundary(), later refactor
  to provide FE vector
*/

namespace aspect
{
  namespace MeshDeformation
  {


    /**
     * Produce graphical output defined in points provided by the user.
     */
    template <int dim, int spacedim = dim>
    class PointDataOut : public dealii::DataOutInterface<0, spacedim>
    {
      public:
        /**
         * Default constructor.
         */
        PointDataOut() = default;

        /**
         * Default destructor.
         */
        ~PointDataOut() = default;


        /**
         * Build the patches for a given set of points and optionally data in each point.
         *
         * @param [in] locations The point locations.
         * @param [in] data A vector of data values for each point.
         * @param [in] data_component_names An optional vector of strings that
         * describe the properties of each datum.
         * @param [in] data_component_interpretations An optional vector that
         * controls if the properties are interpreted as scalars, vectors,
         * or tensors. Has to be of the same length as @p data_component_names.
         */
        void
        build_patches(const std::vector<Point<spacedim>> &locations,
                      const std::vector<std::vector<double>> &data = {},
                      const std::vector<std::string> &data_component_names = {},
                      const std::vector<
                      DataComponentInterpretation::DataComponentInterpretation>
                      &data_component_interpretations_ = {})
        {
          Assert(data_component_names.size() == data_component_interpretations_.size(),
                 ExcMessage(
                   "When calling PointDataOut::build_patches() with data component "
                   "names and interpretations you need to provide as many data component "
                   "names as interpretations. Provide the same name for components that "
                   "belong to a single vector or tensor."));

          Assert(data.size() == 0 || locations.size(),
                 ExcMessage("You need to either provide no data or data for each point."));
          for (const auto &datum : data)
            Assert(datum.size() == data_component_names.size(),
                   ExcMessage("The data provided in each point needs to have the same number "
                              "of components as names were provided."));

          // Prepend the "id" to the data fields provided by the user:
          dataset_names.clear();
          dataset_names.emplace_back("id");
          dataset_names.insert(dataset_names.end(),
                               data_component_names.begin(),
                               data_component_names.end());

          data_component_interpretations.clear();
          data_component_interpretations.emplace_back(
            DataComponentInterpretation::component_is_scalar);
          data_component_interpretations.insert(
            data_component_interpretations.end(),
            data_component_interpretations_.begin(),
            data_component_interpretations_.end());

          const unsigned int n_property_components = data_component_names.size();
          const unsigned int n_data_components     = dataset_names.size();

          patches.resize(locations.size());

          for (unsigned int i = 0; i < locations.size(); ++i)
            {
              patches[i].vertices[0] = locations[i];
              patches[i].patch_index = i;

              // Store id and properties given by the user:
              patches[i].data.reinit(n_data_components, 1);
              patches[i].data(0, 0) = i; // store id
              for (unsigned int property_index = 0; property_index < n_property_components; ++property_index)
                patches[i].data(property_index + 1, 0) = data[i][property_index];
            }
        }

      protected:
        /**
         * Returns the patches previously built by the build_patches() function.
         */
        virtual const std::vector<DataOutBase::Patch<0, spacedim>> &
        get_patches() const override
        {
          return patches;
        }

        /**
         * Virtual function through which the names of data sets are obtained from
         * this class
         */
        virtual std::vector<std::string>
        get_dataset_names() const override
        {
          return dataset_names;
        }


        /**
         * Overload of the respective DataOutInterface::get_nonscalar_data_ranges()
         * function. See there for a more extensive documentation.
         * This function is a reimplementation of the function
         * DataOut_DoFData::get_nonscalar_data_ranges().
         */
        virtual std::vector<
        std::tuple<unsigned int,
            unsigned int,
            std::string,
            DataComponentInterpretation::DataComponentInterpretation>>
            get_nonscalar_data_ranges() const override
        {
          std::vector<
          std::tuple<unsigned int,
              unsigned int,
              std::string,
              DataComponentInterpretation::DataComponentInterpretation>>
              ranges;

          // Make sure the data structures were set up correctly. Since they
          // can only be filled by build_patches() above, they should have
          // been checked already.
          Assert(dataset_names.size() == data_component_interpretations.size(),
                 ExcInternalError());

          // collect the ranges of particle data
          const unsigned int n_output_components =
            data_component_interpretations.size();
          unsigned int output_component = 0;
          for (unsigned int i = 0; i < n_output_components; /* i is updated below */)
            // see what kind of data we have here. note that for the purpose of the
            // current function all we care about is vector data
            switch (data_component_interpretations[i])
              {
                case DataComponentInterpretation::component_is_scalar:
                {
                  // Just move component forward by one
                  ++i;
                  ++output_component;

                  break;
                }
                case DataComponentInterpretation::component_is_part_of_vector:
                {
                  // ensure that there is a continuous number of next space_dim
                  // components that all deal with vectors
                  Assert(
                    i + spacedim <= n_output_components,
                    Exceptions::DataOutImplementation::ExcInvalidVectorDeclaration(
                      i, dataset_names[i]));
                  for (unsigned int dd = 1; dd < spacedim; ++dd)
                    Assert(
                      data_component_interpretations[i + dd] ==
                      DataComponentInterpretation::component_is_part_of_vector,
                      Exceptions::DataOutImplementation::
                      ExcInvalidVectorDeclaration(i, dataset_names[i]));

                  // all seems right, so figure out whether there is a common
                  // name to these components. if not, leave the name empty and
                  // let the output format writer decide what to do here
                  std::string name = dataset_names[i];
                  for (unsigned int dd = 1; dd < spacedim; ++dd)
                    if (name != dataset_names[i + dd])
                      {
                        name = "";
                        break;
                      }

                  // Finally add a corresponding range.
                  //
                  // This sort of logic is also explained in some detail in
                  //   DataOut::build_one_patch().
                  ranges.emplace_back(std::forward_as_tuple(
                                        output_component,
                                        output_component + spacedim - 1,
                                        name,
                                        DataComponentInterpretation::component_is_part_of_vector));

                  // increase the 'component' counter by the appropriate amount,
                  // same for 'i', since we have already dealt with all these
                  // components
                  output_component += spacedim;
                  i += spacedim;

                  break;
                }

                case DataComponentInterpretation::component_is_part_of_tensor:
                {
                  const unsigned int size = spacedim * spacedim;
                  // ensure that there is a continuous number of next
                  // spacedim*spacedim components that all deal with tensors
                  Assert(
                    i + size <= n_output_components,
                    Exceptions::DataOutImplementation::ExcInvalidTensorDeclaration(
                      i, dataset_names[i]));
                  for (unsigned int dd = 1; dd < size; ++dd)
                    Assert(
                      data_component_interpretations[i + dd] ==
                      DataComponentInterpretation::component_is_part_of_tensor,
                      Exceptions::DataOutImplementation::
                      ExcInvalidTensorDeclaration(i, dataset_names[i]));

                  // all seems right, so figure out whether there is a common
                  // name to these components. if not, leave the name empty and
                  // let the output format writer decide what to do here
                  std::string name = dataset_names[i];
                  for (unsigned int dd = 1; dd < size; ++dd)
                    if (name != dataset_names[i + dd])
                      {
                        name = "";
                        break;
                      }

                  // Finally add a corresponding range.
                  ranges.emplace_back(std::forward_as_tuple(
                                        output_component,
                                        output_component + size - 1,
                                        name,
                                        DataComponentInterpretation::component_is_part_of_tensor));

                  // increase the 'component' counter by the appropriate amount,
                  // same for 'i', since we have already dealt with all these
                  // components
                  output_component += size;
                  i += size;
                  break;
                }

                default:
                  Assert(false, ExcNotImplemented());
              }
          return ranges;
        }

      private:
        /**
         * This is a vector of patches that is created each time build_patches() is
         * called. These patches are used in the output routines of the base
         * classes.
         */
        std::vector<DataOutBase::Patch<0, spacedim>> patches;

        /**
         * A vector of field names for all data components stored in patches.
         */
        std::vector<std::string> dataset_names;

        /**
         * A vector that for each of the data components of the
         * current data set indicates whether they are scalar fields, parts of a
         * vector-field, or any of the other supported kinds of data.
         */
        std::vector<DataComponentInterpretation::DataComponentInterpretation>
        data_component_interpretations;
    };




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

      {
        // Visualize the evaluation points and their velocities
        static unsigned int output_no = 0;

        PointDataOut<dim, dim> out;
        // const auto &mapping = this->get_mapping();
        std::vector<Point<dim>> real_evaluation_points(evaluation_points.size());
        std::vector<std::vector<double>> data(evaluation_points.size(), std::vector<double>(dim, 0.0));
        for (unsigned int i=0; i<evaluation_points.size(); ++i)
          {
            real_evaluation_points[i] = evaluation_points[i];  // TODO: use mapping to compute real position
            for (unsigned int c=0; c<dim; ++c)
              data[i][c] = velocities[i][c];
          }

        const std::vector<std::string> data_component_names(dim, "velocity");
        const std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretations(dim, DataComponentInterpretation::component_is_part_of_vector);

        out.build_patches(real_evaluation_points, data, data_component_names, data_component_interpretations);

        out.write_vtu_with_pvtu_record(this->get_output_directory(), "surf_points", output_no, this->get_mpi_communicator(), 4, 0);

        ++output_no;
      }



      // Create the output vector.
      const DoFHandler<dim> &mesh_dof_handler = this->get_mesh_deformation_handler().get_mesh_deformation_dof_handler();
      LinearAlgebra::Vector vector_with_surface_velocities(mesh_dof_handler.locally_owned_dofs(),
                                                           this->get_mpi_communicator());

      // The remote_point_evaluator will gives us the velocities in all evaluation points that are within one of our locally
      // owned cells. The lambda defined below receives a list of points and their velocities for each cell. The coordinates
      // are given in coordinates of the unit cell.
      // For each support point of the velocity DoFHandler, we will set the velocity from the closest evaluation point. We
      // do this by keeping track of 1/distance of the closest evaluation point checked so far. The initial value of 0.0
      // denotes an infinite distance.
      LinearAlgebra::Vector one_over_distance_vec(mesh_dof_handler.locally_owned_dofs(),
                                                  this->get_mpi_communicator());

      const unsigned int dofs_per_cell = mesh_dof_handler.get_fe().dofs_per_cell;

      // Note: We assume that process_and_evaluate() does not call our lambda concurrently, otherwise we would have write
      // conflicts when updating vector_with_surface_velocities and one_over_distance_vec.
      const auto eval_func = [&](const ArrayView< const Tensor<1,dim>> &values,
                                 const typename Utilities::MPI::RemotePointEvaluation<dim>::CellData &cell_data)
      {
        std::vector<types::global_dof_index> cell_dof_indices (dofs_per_cell);
        for (const auto cell_index : cell_data.cell_indices())
          {
            const auto cell_dofs =
              cell_data.get_active_cell_iterator(cell_index)->as_dof_handler_iterator(
                mesh_dof_handler);

            const ArrayView<const Point<dim>> unit_points = cell_data.get_unit_points(cell_index);
            const auto local_values = cell_data.get_data_view(cell_index, values);

            cell_dofs->get_dof_indices(cell_dof_indices);

            // Note: This search is a nested loop with the inner part executed #evaluation_point_in_this_cell * #dofs_per_cell
            // times. We could precompute this information as the point locations and do not change (outside of mesh refinement).
            const std::vector< Point< dim >> &support_points = mesh_dof_handler.get_fe().get_unit_support_points();
            for (unsigned int i=0; i<unit_points.size(); ++i)
              {
                for (unsigned int j=0; j<support_points.size(); ++j)
                  {
                    const double one_over_distance = 1.0/(unit_points[i].distance(support_points[j])+1e-10);
                    if (one_over_distance > one_over_distance_vec(cell_dof_indices[j]))
                      {
                        // The point i is closer to support point j than anything we have seen so far. Keep track
                        // of the distance and write the correct velocity component into the result:
                        one_over_distance_vec(cell_dof_indices[j]) = one_over_distance;
                        const unsigned int component = mesh_dof_handler.get_fe().system_to_component_index(j).first;
                        vector_with_surface_velocities(cell_dof_indices[j]) = local_values[i][component];
                      }
                  }
              }
          }
      };

      this->remote_point_evaluator->template process_and_evaluate<Tensor<1,dim>,1>(velocities, eval_func, /*sort_data*/ true);
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
