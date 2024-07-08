#include <aspect/simulator.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/mesh_deformation/interface.h>

namespace aspect
{
  template <int dim>
  void
  Simulator<dim>::apply_BP_limiter_to_dg_solutions(const AdvectionField &advection_field)
  {
    // TODO: Modify to more robust method
    // Skip if this composition field is being set from the volume_of_fluid handler
    if (!advection_field.is_temperature() &&
        parameters.volume_of_fluid_tracking_enabled)
      if (volume_of_fluid_handler->field_index_for_name(introspection.name_for_compositional_index(advection_field.compositional_variable))
          != volume_of_fluid_handler->get_n_fields())
        return;

    /*
     * First setup the quadrature points which are used to find the maximum and minimum solution values at those points.
     * A quadrature formula that combines all quadrature points constructed as all tensor products of
     * 1) one dimensional Gauss points; 2) one dimensional Gauss-Lobatto points.
     * We require that the Gauss-Lobatto points (2) appear in only one direction.
     * Therefore, possible combination
     * in 2d: the combinations are 21, 12
     * in 3d: the combinations are 211, 121, 112
     */
    const QGauss<1> quadrature_formula_1 (advection_field.polynomial_degree(introspection)+1);
    const QGaussLobatto<1> quadrature_formula_2 (advection_field.polynomial_degree(introspection)+1);

    const unsigned int n_q_points_1 = quadrature_formula_1.size();
    const unsigned int n_q_points_2 = quadrature_formula_2.size();
    const unsigned int n_q_points   = dim * n_q_points_2 * static_cast<unsigned int>(std::pow(n_q_points_1, dim-1));

    std::vector<Point <dim>> quadrature_points;
    quadrature_points.reserve(n_q_points);

    switch (dim)
      {
        case 2:
        {
          // append quadrature points combination 12
          for (unsigned int i=0; i < n_q_points_1 ; ++i)
            {
              const double  x = quadrature_formula_1.point(i)(0);
              for (unsigned int j=0; j < n_q_points_2 ; ++j)
                {
                  const double y = quadrature_formula_2.point(j)(0);
                  quadrature_points.push_back(Point<dim>(x,y));
                }
            }
          // append quadrature points combination 21
          for (unsigned int i=0; i < n_q_points_2 ; ++i)
            {
              const double  x = quadrature_formula_2.point(i)(0);
              for (unsigned int j=0; j < n_q_points_1 ; ++j)
                {
                  const double y = quadrature_formula_1.point(j)(0);
                  quadrature_points.push_back(Point<dim>(x,y));
                }
            }
          break;
        }

        case 3:
        {
          // append quadrature points combination 121
          for ( unsigned int i=0; i < n_q_points_1 ; ++i)
            {
              const double x = quadrature_formula_1.point(i)(0);
              for ( unsigned int j=0; j < n_q_points_2 ; ++j)
                {
                  const double y = quadrature_formula_2.point(j)(0);
                  for ( unsigned int k=0; k < n_q_points_1 ; ++k)
                    {
                      const double z = quadrature_formula_1.point(k)(0);
                      quadrature_points.push_back(Point<dim>(x,y,z));
                    }
                }
            }
          // append quadrature points combination 112
          for (unsigned int i=0; i < n_q_points_1 ; ++i)
            {
              const double x = quadrature_formula_1.point(i)(0);
              for (unsigned int j=0; j < n_q_points_1 ; ++j)
                {
                  const double y = quadrature_formula_1.point(j)(0);
                  for (unsigned int k=0; k < n_q_points_2 ; ++k)
                    {
                      const double z = quadrature_formula_2.point(k)(0);
                      quadrature_points.push_back(Point<dim>(x,y,z));
                    }
                }
            }
          // append quadrature points combination 211
          for (unsigned int i=0; i < n_q_points_2 ; ++i)
            {
              const double x = quadrature_formula_2.point(i)(0);
              for ( unsigned int j=0; j < n_q_points_1 ; ++j)
                {
                  const double y = quadrature_formula_1.point(j)(0);
                  for ( unsigned int k=0; k < n_q_points_1 ; ++k)
                    {
                      const double z = quadrature_formula_1.point(k)(0);
                      quadrature_points.push_back(Point<dim>(x,y,z));
                    }
                }
            }
          break;
        }

        default:
          Assert (false, ExcNotImplemented());
      }

    Assert (quadrature_points.size() == n_q_points, ExcInternalError());
    const Quadrature<dim> quadrature_formula(quadrature_points);

    // Quadrature rules only used for the numerical integration for better accuracy
    const Quadrature<dim> &quadrature_formula_0
      = (advection_field.is_temperature() ?
         introspection.quadratures.temperature :
         introspection.quadratures.compositional_fields);
    const unsigned int n_q_points_0 = quadrature_formula_0.size();

    // fe values for points evaluation
    FEValues<dim> fe_values (*mapping,
                             finite_element,
                             quadrature_formula,
                             update_values   |
                             update_quadrature_points);
    std::vector<double> values (n_q_points);
    // fe values for numerical integration, with a number of quadrature points
    // that is equal to 1/dim times the number of total points above
    FEValues<dim> fe_values_0 (*mapping,
                               finite_element,
                               quadrature_formula_0,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
    std::vector<double> values_0 (n_q_points_0);

    const FEValuesExtractors::Scalar field
      = (advection_field.is_temperature()
         ?
         introspection.extractors.temperature
         :
         introspection.extractors.compositional_fields[advection_field.compositional_variable]
        );

    const double max_solution_exact_global = (advection_field.is_temperature()
                                              ?
                                              parameters.global_temperature_max_preset
                                              :
                                              parameters.global_composition_max_preset[advection_field.compositional_variable]
                                             );
    const double min_solution_exact_global = (advection_field.is_temperature()
                                              ?
                                              parameters.global_temperature_min_preset
                                              :
                                              parameters.global_composition_min_preset[advection_field.compositional_variable]
                                             );

    LinearAlgebra::BlockVector distributed_solution (introspection.index_sets.system_partitioning,
                                                     mpi_communicator);
    const unsigned int block_idx = advection_field.block_index(introspection);
    distributed_solution.block(block_idx) = solution.block(block_idx);

    std::vector<types::global_dof_index> local_dof_indices (finite_element.dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (local_dof_indices);
          // used to find the maximum, minimum
          fe_values.reinit (cell);
          fe_values[field].get_function_values(solution, values);
          // used for the numerical integration
          fe_values_0.reinit (cell);
          fe_values_0[field].get_function_values(solution, values_0);

          // Find the local max and local min
          const double min_solution_local = *std::min_element (values.begin(), values.end());
          const double max_solution_local = *std::max_element (values.begin(), values.end());
          // Find the trouble cell
          if (min_solution_local < min_solution_exact_global
              || max_solution_local > max_solution_exact_global)
            {
              // Compute the cell area and cell solution average
              double local_area = 0.0;
              double local_solution_average = 0.0;
              for (unsigned int q = 0; q < n_q_points_0; ++q)
                {
                  local_area += fe_values_0.JxW(q);
                  local_solution_average += values_0[q]*fe_values_0.JxW(q);
                }
              local_solution_average /= local_area;

              /*
               * Define theta: a scaling constant used to correct the old solution by the formula
               *   new_value = theta * (old_value-old_solution_cell_average)+old_solution_cell_average
               * where theta \in [0,1] defined as below.
               * After the correction, the new solution does not exceed the user-given
               * exact global maximum/minimum values. Meanwhile, the new solution's cell average
               * equals to the old solution's cell average.
               */
              double theta = 1.0;
              if (std::abs(max_solution_local-local_solution_average) > std::numeric_limits<double>::min())
                {
                  theta = std::min(theta, std::abs((max_solution_exact_global-local_solution_average)
                                                   / (max_solution_local-local_solution_average)));
                }
              if (std::abs(min_solution_local-local_solution_average) > std::numeric_limits<double>::min())
                {
                  theta = std::min(theta, std::abs((min_solution_exact_global-local_solution_average)
                                                   / (min_solution_local-local_solution_average)));
                }

              /* Modify the advection degrees of freedom of the numerical solution.
               * Note that we are using DG elements, so every DoF on a locally owned cell is locally owned;
               * this means that we do not need to check whether the 'distributed_solution' vector actually
               * stores the element we read from/write to here.
               */
              for (unsigned int j = 0;
                   j < finite_element.base_element(advection_field.base_element(introspection)).dofs_per_cell;
                   ++j)
                {
                  const unsigned int support_point_index = finite_element.component_to_system_index(
                                                             (advection_field.is_temperature()
                                                              ?
                                                              introspection.component_indices.temperature
                                                              :
                                                              introspection.component_indices.compositional_fields[advection_field.compositional_variable]
                                                             ),
                                                             /*dof index within component=*/ j);
                  const double solution_value = solution(local_dof_indices[support_point_index]);
                  const double limited_solution_value = theta * (solution_value-local_solution_average) + local_solution_average;
                  distributed_solution(local_dof_indices[support_point_index]) = limited_solution_value;
                }
            }
        }

    distributed_solution.compress(VectorOperation::insert);
    // now get back to the original vector
    solution.block(block_idx) = distributed_solution.block(block_idx);
  }



  namespace internal
  {
    struct KXRCFCellData
    {
      double cell_norm;
      double inflow_face_jump;
      double inflow_face_area;

      KXRCFCellData()
        : cell_norm(0.0)
        , inflow_face_jump(0.0)
        , inflow_face_area(0.0)
      {}
    };
  }


  template <int dim>
  template <typename T>
  void
  Simulator<dim>::compute_KXRCF_indicators(Vector<T>            &KXRCF_indicators,
                                           const AdvectionField &advection_field) const
  {
    // stuff for computing the integrals in cells and cell faces
    QGauss<dim>   quadrature(2);
    QGauss<dim-1> face_quadrature(2);

    FEValues<dim> fe_values(*mapping,
                            finite_element,
                            quadrature,
                            update_values);

    FEFaceValues<dim> fe_face_values(*mapping,
                                     finite_element,
                                     face_quadrature,
                                     update_values |
                                     update_normal_vectors |
                                     update_JxW_values);

    FESubfaceValues<dim> fe_subface_values(*mapping,
                                           finite_element,
                                           face_quadrature,
                                           update_values |
                                           update_normal_vectors |
                                           update_JxW_values);

    FEFaceValues<dim> neighbor_fe_face_values(*mapping,
                                              finite_element,
                                              face_quadrature,
                                              update_values |
                                              update_JxW_values);

    const FEValuesExtractors::Scalar &field_extractor = advection_field.scalar_extractor(introspection);

    const unsigned int n_q_points      = fe_values.n_quadrature_points;
    const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

    std::vector<double> field_values(n_q_points);
    std::vector<double> face_field_values(n_face_q_points);
    std::vector<double> neighbor_face_field_values(n_face_q_points);

    std::vector<Tensor<1,dim>> face_velocity_values(n_face_q_points);
    std::vector<Tensor<1,dim>> face_mesh_velocity_values(n_face_q_points);

    // KXRCF indicator requires the computation of (1) cell maximum norm;
    // (2) integration of the jump over inflow faces of each field. To
    // avoid doing integration twice on each cell face, we set up a map
    // to store the cell data and then loop over the cells using the same
    // strategy as in Assemblers::AdvectionSystemInteriorFace::execute().
    std::map<unsigned int, internal::KXRCFCellData> data;
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        data.insert(std::make_pair(cell->active_cell_index(),
                                   internal::KXRCFCellData()));

    // Now loop over locally owned cells and compute the maximum norm
    // and the jump over inflow faces.
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          auto cell_data = data.find(cell->active_cell_index());

          // Compute the maximum norm of field values in this cell.
          fe_values.reinit(cell);
          fe_values[field_extractor].get_function_values(solution, field_values);

          double cell_norm = 0.0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            cell_norm = std::max(cell_norm, std::abs(field_values[q]));

          cell_data->second.cell_norm = cell_norm;

          // Compute the jump of field values over inflow faces.
          for (const unsigned int face_no : cell->face_indices())
            {
              if (!cell->at_boundary(face_no) || cell->has_periodic_neighbor(face_no))
                {
                  const typename DoFHandler<dim>::cell_iterator
                  neighbor = cell->neighbor_or_periodic_neighbor(face_no);

                  const bool cell_has_periodic_neighbor = cell->has_periodic_neighbor(face_no);

                  // 'neighbor' defined above is NOT active_cell_iterator, this includes cells
                  // that are refined
                  if (!neighbor->has_children())
                    {
                      if (neighbor->level() == cell->level() &&
                          neighbor->is_active() &&
                          (((neighbor->is_locally_owned()) && (cell->index() < neighbor->index()))
                           ||
                           ((!neighbor->is_locally_owned()) && (cell->subdomain_id() < neighbor->subdomain_id()))))
                        {
                          // cell and neighbor are equal-sized, and cell has been chosen to
                          // assemble this face, so calculate from cell
                          const unsigned int neighbor2 =
                            (cell_has_periodic_neighbor
                             ?
                             cell->periodic_neighbor_of_periodic_neighbor(face_no)
                             :
                             cell->neighbor_of_neighbor(face_no));

                          auto neighbor_data = data.find(neighbor->active_cell_index());

                          // Set up face values.
                          fe_face_values.reinit(cell, face_no);
                          fe_face_values[field_extractor].get_function_values(solution, face_field_values);

                          fe_face_values[introspection.extractors.velocities].get_function_values(
                            solution, face_velocity_values);
                          if (parameters.mesh_deformation_enabled)
                            fe_face_values[introspection.extractors.velocities].get_function_values(
                              mesh_deformation->mesh_velocity, face_mesh_velocity_values);

                          // Set up neighbor face values.
                          neighbor_fe_face_values.reinit(neighbor, neighbor2);
                          neighbor_fe_face_values[field_extractor].get_function_values(solution, neighbor_face_field_values);

                          const double face_area = cell->face(face_no)->measure();

                          double v_n = 0.0;
                          double face_jump = 0.0;
                          for (unsigned int q = 0; q < n_face_q_points; ++q)
                            {
                              Tensor<1,dim> current_u = face_velocity_values[q];
                              if (parameters.mesh_deformation_enabled)
                                current_u -= face_mesh_velocity_values[q];

                              v_n += (current_u * fe_face_values.normal_vector(q)) * fe_face_values.JxW(q);

                              face_jump += (face_field_values[q] - neighbor_face_field_values[q]) * fe_face_values.JxW(q);
                            }

                          if (v_n < 0.0) // inflow
                            {
                              cell_data->second.inflow_face_jump += face_jump;
                              cell_data->second.inflow_face_area += face_area;
                            }
                          else // outflow
                            {
                              if (neighbor_data != data.end())
                                {
                                  neighbor_data->second.inflow_face_jump -= face_jump;
                                  neighbor_data->second.inflow_face_area += face_area;
                                }
                            }

                        }
                      else
                        {
                          // Neighbor is taking responsibility for assembly of this face, because
                          // either (1) neighbor is coarser, or
                          //        (2) neighbor is equally-sized and
                          //           (a) neighbor is on a different subdomain, with lower subdomain_id(), or
                          //           (b) neighbor is on the same subdomain and has lower index().
                        }
                    }
                  else
                    {
                      // Neighbor has children, so always assemble from here.
                      const unsigned int neighbor2 =
                        (cell_has_periodic_neighbor
                         ?
                         cell->periodic_neighbor_face_no(face_no)
                         :
                         cell->neighbor_face_no(face_no));

                      // Loop over subfaces. We know that the neighbor is finer, so we could loop over the subfaces of the current
                      // face. but if we are at a periodic boundary, then the face of the current cell has no children, so instead use
                      // the children of the periodic neighbor's corresponding face since we know that the letter does indeed have
                      // children (because we know that the neighbor is refined).
                      typename DoFHandler<dim>::face_iterator neighbor_face = neighbor->face(neighbor2);
                      for (unsigned int subface_no = 0; subface_no < neighbor_face->n_children(); ++subface_no)
                        {
                          const typename DoFHandler<dim>::active_cell_iterator neighbor_child =
                            (cell_has_periodic_neighbor
                             ?
                             cell->periodic_neighbor_child_on_subface(face_no, subface_no)
                             :
                             cell->neighbor_child_on_subface(face_no, subface_no));

                          auto neighbor_data = data.find(neighbor_child->active_cell_index());

                          // Set up subface values.
                          fe_subface_values.reinit(cell, face_no, subface_no);
                          fe_subface_values[field_extractor].get_function_values(solution, face_field_values);

                          fe_subface_values[introspection.extractors.velocities].get_function_values(
                            solution, face_velocity_values);
                          if (parameters.mesh_deformation_enabled)
                            fe_subface_values[introspection.extractors.velocities].get_function_values(
                              mesh_deformation->mesh_velocity, face_mesh_velocity_values);

                          neighbor_fe_face_values.reinit(neighbor_child, neighbor2);
                          neighbor_fe_face_values[field_extractor].get_function_values(solution, neighbor_face_field_values);

                          const double face_area = neighbor_child->face(neighbor2)->measure();

                          double v_n = 0.0;
                          double face_jump = 0.0;
                          for (unsigned int q = 0; q < n_face_q_points; ++q)
                            {
                              Tensor<1,dim> current_u = face_velocity_values[q];
                              if (parameters.mesh_deformation_enabled)
                                current_u -= face_mesh_velocity_values[q];

                              v_n += (current_u * fe_subface_values.normal_vector(q)) * fe_subface_values.JxW(q);

                              face_jump += (face_field_values[q] - neighbor_face_field_values[q]) * fe_subface_values.JxW(q);
                            }

                          if (v_n < 0.0) // inflow
                            {
                              cell_data->second.inflow_face_jump += face_jump;
                              cell_data->second.inflow_face_area += face_area;
                            }
                          else // outflow
                            {
                              if (neighbor_data != data.end())
                                {
                                  neighbor_data->second.inflow_face_jump -= face_jump;
                                  neighbor_data->second.inflow_face_area += face_area;
                                }
                            }
                        }
                    }
                }
              else
                {
                  // The current face is a boundary face.
                }
            }
        }

    // Finally, compute the KXRCF indicator and pick the largest one as
    // troubled cell indicator.
    const unsigned int degree = advection_field.polynomial_degree(introspection);
    const double power = (degree + 1.0) / 2.0;
    const double scaling_factor = 1. / global_Omega_diameter;

    const unsigned int block_idx = advection_field.block_index(introspection);
    const double epsilon = std::max(1e-20,
                                    1e-3 * (solution.block(block_idx).max() - solution.block(block_idx).min()));

    KXRCF_indicators.reinit(triangulation.n_active_cells());
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          // The characteristic mesh spacing is calculated by d/D, where
          // d and D are the diameters of the present cell and the geometry
          // model, respectively.
          const double h_pow = std::pow(cell->diameter() * scaling_factor, power);
          const auto cell_data = data.find(cell->active_cell_index());
          if (cell_data->second.inflow_face_area == 0.0)
            KXRCF_indicators[cell->active_cell_index()] = 0.0;
          else
            KXRCF_indicators[cell->active_cell_index()] =
              std::abs(cell_data->second.inflow_face_jump) /
              (h_pow * cell_data->second.inflow_face_area * (cell_data->second.cell_norm + epsilon));
        }
  }


  template <int dim>
  void
  Simulator<dim>::apply_WENO_limiter_to_dg_solutions(const AdvectionField &advection_field)
  {
    // Skip if this composition field is being set from the volume_of_fluid handler
    if (!advection_field.is_temperature() &&
        parameters.volume_of_fluid_tracking_enabled)
      if (volume_of_fluid_handler->field_index_for_name(introspection.name_for_compositional_index(advection_field.compositional_variable))
          != volume_of_fluid_handler->get_n_fields())
        return;

    const unsigned int degree = advection_field.polynomial_degree(introspection);
    AssertThrow(degree < 3,
                ExcMessage("Currently, WENO limiter has only been implemented for "
                           "advection fields with polynomial degree 1 or 2."));

    // We need two objects of FEValues in the following computation:
    // fe_values: computes the cell average and smoothness indicator of the target cell,
    //   and provides the coordinates and JxWs of the quadrature points of the target
    //   cell when computing the smoothness indicators of neighbor cells;
    // neighbor_fe_values: computes the cell averages and field derivatives of
    //   neighbor cells at quadrature points of the target cell. With the derivatives
    //   at hand, we can use the JxWs provided by fe_values to compute the smoothness
    //   indicators of neighbor cells.
    // The coordinates of quadrature points of neighbor_fe_values change from cell to
    // cell since the mesh may not be Cartesian, hence neighbor_fe_values has to be
    // reconstructed in each cell.
    Quadrature<dim> quadrature(advection_field.is_temperature() ?
                               introspection.quadratures.temperature :
                               introspection.quadratures.compositional_fields);

    UpdateFlags update_flags = update_values | update_gradients |
                               update_quadrature_points | update_JxW_values |
                               (degree > 1 ? update_hessians : update_default);

    UpdateFlags neighbor_update_flags = update_values | update_gradients |
                                        (degree > 1 ? update_hessians : update_default);

    FEValues<dim> fe_values(*mapping, finite_element, quadrature, update_flags);

    const unsigned int n_q_points          = fe_values.n_quadrature_points;
    const unsigned int field_dofs_per_cell = finite_element.base_element(advection_field.base_element(introspection)).dofs_per_cell;
    const unsigned int field_component     = advection_field.component_index(introspection);

    const FEValuesExtractors::Scalar &field_extractor = advection_field.scalar_extractor(introspection);

    std::vector<double>        field_values(n_q_points);
    std::vector<Tensor<1,dim>> field_gradients(n_q_points);
    std::vector<Tensor<2,dim>> field_hessians(n_q_points);

    std::vector<Point<dim>> neighbor_quadrature_points(n_q_points);

    std::vector<double> reconstructed_values(n_q_points);

    std::vector<double> cell_averages;
    std::vector<double> smoothness_indicators;
    std::vector<double> linear_weights;
    std::vector<double> nonlinear_weights;
    std::vector<std::vector<double>> field_values_in_target_cell;

    // stuff for projecting the reconstructed polynomial onto grid nodes
    FullMatrix<double> cell_matrix(field_dofs_per_cell, field_dofs_per_cell);
    Vector<double> cell_rhs(field_dofs_per_cell);
    Vector<double> cell_solution(field_dofs_per_cell);
    std::vector<double> phi_field(field_dofs_per_cell);
    std::vector<types::global_dof_index> cell_dof_indices(finite_element.dofs_per_cell);

    // Detect troubled cells with KXRCF indicator.
    Vector<double> KXRCF_indicators;
    compute_KXRCF_indicators(KXRCF_indicators, advection_field);

    TrilinosWrappers::MPI::BlockVector distributed_solution(introspection.index_sets.system_partitioning,
                                                            mpi_communicator);

    const unsigned int block_idx = advection_field.block_index(introspection);
    distributed_solution.block(block_idx) = solution.block(block_idx);

    const double KXRCF_indicator_threshold = (advection_field.is_temperature() ?
                                              parameters.temperature_KXRCF_indicator_threshold :
                                              parameters.composition_KXRCF_indicator_threshold[advection_field.compositional_variable]);
    for (const auto &cell : dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          // If the KXRCF indicator of the present cell is lower than the threshold,
          // then we do not have to do WENO reconstruction for the cell.
          if (KXRCF_indicators[cell->active_cell_index()] < KXRCF_indicator_threshold)
            continue;

          cell_averages.clear();
          smoothness_indicators.clear();
          field_values_in_target_cell.clear();
          linear_weights.clear();

          // Compute the field values and derivatives at quadrature points.
          fe_values.reinit(cell);
          fe_values[field_extractor].get_function_values(solution, field_values);
          fe_values[field_extractor].get_function_gradients(solution, field_gradients);
          if (degree > 1)
            fe_values[field_extractor].get_function_hessians(solution, field_hessians);

          // Compute the cell average and the smoothness indicator.
          double cell_measure                   = 0.0;
          double field_integral                 = 0.0;
          double field_gradient_square_integral = 0.0;
          double field_hessian_square_integral  = 0.0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              cell_measure += fe_values.JxW(q);
              field_integral += field_values[q] * fe_values.JxW(q);

              field_gradient_square_integral += field_gradients[q].norm_square() * fe_values.JxW(q);
              if (degree > 1)
                field_hessian_square_integral += field_hessians[q].norm_square() * fe_values.JxW(q);
            }

          double beta = field_gradient_square_integral;
          if (degree > 1)
            beta += field_hessian_square_integral * cell_measure;

          cell_averages.push_back(field_integral / cell_measure);
          smoothness_indicators.push_back(beta);
          field_values_in_target_cell.push_back(field_values);
          // Linear weight of the target cell cannot be determined at present time
          // (because we don't know how many neighbors the target cell has). Fill in
          // the position with an invalid value.
          linear_weights.push_back(numbers::signaling_nan<double>());

          // Loop over the neighbor cells.
          for (const unsigned int face_no : cell->face_indices())
            {
              if (cell->face(face_no)->at_boundary())
                continue;

              const typename DoFHandler<dim>::cell_iterator
              neighbor = cell->neighbor_or_periodic_neighbor(face_no);

              // 'neighbor' defined above is NOT active_cell_iterator, this includes cells
              // that are refined
              if (!neighbor->has_children())
                {
                  // Compute the field values at quadrature points of the neighbor.

                  // Compute the field derivatives at quadrature points of the target cell.
                  // To do so, we need to transform the quadrature points of the target cell
                  // to the 'unit cell' defined by the neighbor.
                  mapping->transform_points_real_to_unit_cell(neighbor,
                                                              fe_values.get_quadrature_points(),
                                                              neighbor_quadrature_points);

                  Quadrature<dim> neighbor_quadrature(neighbor_quadrature_points,
                                                      quadrature.get_weights());

                  FEValues<dim> neighbor_fe_values(*mapping,
                                                   finite_element,
                                                   neighbor_quadrature,
                                                   neighbor_update_flags);

                  neighbor_fe_values.reinit(neighbor);
                  neighbor_fe_values[field_extractor].get_function_values(solution, field_values);
                  neighbor_fe_values[field_extractor].get_function_gradients(solution, field_gradients);
                  if (degree > 1)
                    neighbor_fe_values[field_extractor].get_function_hessians(solution, field_hessians);

                  // Compute the cell average and the smoothness indicator.
                  double field_integral                 = 0.0;
                  double field_gradient_square_integral = 0.0;
                  double field_hessian_square_integral  = 0.0;
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    {
                      field_integral += field_values[q] * fe_values.JxW(q);
                      field_gradient_square_integral += field_gradients[q].norm_square() * fe_values.JxW(q);
                      if (degree > 1)
                        field_hessian_square_integral += field_hessians[q].norm_square() * fe_values.JxW(q);
                    }

                  double beta = field_gradient_square_integral;
                  if (degree > 1)
                    beta += field_hessian_square_integral * cell_measure;

                  cell_averages.push_back(field_integral / cell_measure);
                  smoothness_indicators.push_back(beta);
                  field_values_in_target_cell.push_back(field_values);
                  linear_weights.push_back(parameters.WENO_linear_weight);
                }
              else
                {
                  // Neighbor has children. We use a half of the children that are adjacent to
                  // the target cell (i.e., children that share a face with the target cell) to
                  // compute the cell average and smoothness indicator.
                  const unsigned int neighbor_face_no =
                    (cell->has_periodic_neighbor(face_no)
                     ?
                     cell->periodic_neighbor_face_no(face_no)
                     :
                     cell->neighbor_face_no(face_no));

                  typename DoFHandler<dim>::face_iterator neighbor_face = neighbor->face(neighbor_face_no);
                  for (unsigned int subface_no = 0; subface_no < neighbor_face->n_children(); ++subface_no)
                    {
                      const unsigned int child_no =
                        GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>::isotropic_refinement,
                                                              neighbor_face_no,
                                                              subface_no,
                                                              neighbor->face_orientation(face_no),
                                                              neighbor->face_flip(face_no),
                                                              neighbor->face_rotation(face_no));

                      typename DoFHandler<dim>::active_cell_iterator
                      neighbor_child = neighbor->child(child_no);

                      // Compute the field derivatives at quadrature points of the target cell.
                      // To do so, we need to transform the quadrature points of the target cell
                      // to the 'unit cell' defined by the neighbor.
                      mapping->transform_points_real_to_unit_cell(neighbor_child,
                                                                  fe_values.get_quadrature_points(),
                                                                  neighbor_quadrature_points);

                      Quadrature<dim> neighbor_quadrature(neighbor_quadrature_points,
                                                          quadrature.get_weights());

                      FEValues<dim> neighbor_fe_values(*mapping,
                                                       finite_element,
                                                       neighbor_quadrature,
                                                       neighbor_update_flags);

                      neighbor_fe_values.reinit(neighbor_child);
                      neighbor_fe_values[field_extractor].get_function_values(solution, field_values);
                      neighbor_fe_values[field_extractor].get_function_gradients(solution, field_gradients);
                      if (degree > 1)
                        neighbor_fe_values[field_extractor].get_function_hessians(solution, field_hessians);

                      // Compute the cell average and the smoothness indicator.
                      double field_integral                 = 0.0;
                      double field_gradient_square_integral = 0.0;
                      double field_hessian_square_integral  = 0.0;
                      for (unsigned int q = 0; q < n_q_points; ++q)
                        {
                          field_integral += field_values[q] * fe_values.JxW(q);
                          field_gradient_square_integral += field_gradients[q].norm_square() * fe_values.JxW(q);
                          if (degree > 1)
                            field_hessian_square_integral += field_hessians[q].norm_square() * fe_values.JxW(q);
                        }

                      double beta = field_gradient_square_integral;
                      if (degree > 1)
                        beta += field_hessian_square_integral * cell_measure;

                      cell_averages.push_back(field_integral / cell_measure);
                      smoothness_indicators.push_back(beta);
                      field_values_in_target_cell.push_back(field_values);
                      linear_weights.push_back(parameters.WENO_linear_weight / neighbor_face->n_children());
                    }
                }
            }

          // Compute the linear weights of the target cell.
          linear_weights[0] = 1.0;
          for (unsigned int i = 1; i < linear_weights.size(); ++i)
            linear_weights[0] -= linear_weights[i];

          // Compute the nonlinear weights.
          const double epsilon = 1e-6 * std::accumulate(smoothness_indicators.begin(),
                                                        smoothness_indicators.end(),
                                                        0.0);
          nonlinear_weights.clear();
          for (unsigned int i = 0; i < linear_weights.size(); ++i)
            nonlinear_weights.push_back(linear_weights[i] /
                                        Utilities::fixed_power<2,double>(smoothness_indicators[i] + epsilon));

          // Normalize the nonlinear weights.
          const double sum = std::accumulate(nonlinear_weights.begin(),
                                             nonlinear_weights.end(),
                                             0.0);

          for (unsigned int i = 0; i < nonlinear_weights.size(); ++i)
            nonlinear_weights[i] /= sum;

          // Compute the values of the reconstructed polynomial at quadrature points.
          std::fill(reconstructed_values.begin(), reconstructed_values.end(), 0.0);
          for (unsigned int q = 0; q < n_q_points; ++q)
            for (unsigned int i = 0; i < nonlinear_weights.size(); ++i)
              reconstructed_values[q] += (field_values_in_target_cell[i][q] - cell_averages[i] + cell_averages[0]) * nonlinear_weights[i];

          // Project the reconstructed polynomial onto grid nodes.
          cell_matrix = 0;
          cell_rhs = 0;
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
                {
                  if (finite_element.system_to_component_index(i).first == field_component)
                    {
                      phi_field[i_field] = fe_values[field_extractor].value(i, q);
                      ++i_field;
                    }
                  ++i;
                }

              for (unsigned int i = 0; i < field_dofs_per_cell; ++i)
                {
                  cell_rhs(i) += phi_field[i] * reconstructed_values[q] * fe_values.JxW(q);

                  for (unsigned int j = 0; j < field_dofs_per_cell; ++j)
                    cell_matrix(i, j) += phi_field[i] * phi_field[j] * fe_values.JxW(q);
                }
            }

          cell_matrix.gauss_jordan();
          cell_matrix.vmult(cell_solution, cell_rhs);

          // Copy the results to the distributed solution vector.
          cell->get_dof_indices(cell_dof_indices);
          for (unsigned int i = 0, i_field = 0; i_field < field_dofs_per_cell; /*increment at end of loop*/)
            {
              if (finite_element.system_to_component_index(i).first == field_component)
                {
                  distributed_solution[cell_dof_indices[i]] = cell_solution(i_field);
                  ++i_field;
                }
              ++i;
            }
        }

    // Update the solution vector.
    distributed_solution.compress(VectorOperation::insert);
    solution.block(block_idx) = distributed_solution.block(block_idx);
  }
}

// explicit instantiations
namespace aspect
{

#define INSTANTIATE(dim) \
  template void Simulator<dim>::apply_BP_limiter_to_dg_solutions(const AdvectionField &advection_field); \
  template void Simulator<dim>::compute_KXRCF_indicators(Vector<double> &, \
                                                         const AdvectionField &) const; \
  template void Simulator<dim>::compute_KXRCF_indicators(Vector<float> &, \
                                                         const AdvectionField &) const; \
  template void Simulator<dim>::apply_WENO_limiter_to_dg_solutions(const AdvectionField &advection_field);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
