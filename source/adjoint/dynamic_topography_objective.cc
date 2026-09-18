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

#include <aspect/adjoint/dynamic_topography_objective.h>
#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    DynamicTopographyObjective<dim>::DynamicTopographyObjective ()
      :
      target_topography_model ("none"),
      target_topography_function (1),
      target_coordinate_system (Utilities::Coordinates::spherical)
    {}

    template <int dim>
    double
    DynamicTopographyObjective<dim>::target_topography (const Point<dim> &position) const
    {
      if (target_topography_model == "function")
        {
          const Utilities::NaturalCoordinate<dim> point =
            this->get_geometry_model().cartesian_to_other_coordinates(position,
                                                                      target_coordinate_system);
          return target_topography_function.value(
                   Utilities::convert_array_to_point<dim>(point.get_coordinates()));
        }

      AssertThrow (false,
                   ExcMessage("Dynamic topography objective requires a target topography model. "
                              "Set Adjoint/Dynamic topography/Target topography model."));
      return 0.0;
    }



    template <int dim>
    std::set<std::string>
    DynamicTopographyObjective<dim>::required_forward_postprocessors() const
    {
      return {"dynamic topography"};
    }



    template <int dim>
    double
    DynamicTopographyObjective<dim>::evaluate() const
    {
      const Postprocess::DynamicTopography<dim> &dynamic_topography =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

      AssertThrow(dynamic_topography.topography_vector().size() != 0,
                  ExcMessage("DynamicTopographyObjective requires the dynamic topography "
                             "postprocessor to have run before objective evaluation."));

      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula_face,
                                       update_values |
                                       update_quadrature_points |
                                       update_JxW_values);

      double local_objective = 0.0;
      std::vector<double> topo_values(quadrature_formula_face.size());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_face_values.reinit(cell, f);
                fe_face_values[this->introspection().extractors.temperature].get_function_values(dynamic_topography.topography_vector(),
                    topo_values);

                for (unsigned int q = 0; q < quadrature_formula_face.size(); ++q)
                  {
                    const double residual =
                      topo_values[q] - target_topography(fe_face_values.quadrature_point(q));
                    local_objective += 0.5 * residual * residual * fe_face_values.JxW(q);
                  }
              }

      return Utilities::MPI::sum(local_objective, this->get_mpi_communicator());
    }



    template <int dim>
    void
    DynamicTopographyObjective<dim>::assemble_adjoint_rhs(LinearAlgebra::BlockVector &rhs) const
    {
      const Postprocess::DynamicTopography<dim> &dynamic_topography =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

      AssertThrow(dynamic_topography.topography_vector().size() != 0,
                  ExcMessage("DynamicTopographyObjective requires the dynamic topography "
                             "postprocessor to have run before assembling the adjoint RHS."));

      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> quadrature_formula_face(quadrature_degree);
      const QGaussLobatto<dim-1> cbf_face_quadrature_formula(quadrature_degree);
      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula_face,
                                       update_values |
                                       update_gradients |
                                       update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

      FEFaceValues<dim> fe_cbf_face_values(this->get_mapping(),
                                           this->get_fe(),
                                           cbf_face_quadrature_formula,
                                           update_values |
                                           update_JxW_values);

      const std::vector<Point<dim-1>> face_support_points =
        this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_face_support_points();
      const Quadrature<dim-1> support_quadrature(face_support_points);
      FEFaceValues<dim> fe_support_values(this->get_mapping(),
                                          this->get_fe(),
                                          support_quadrature,
                                          update_values |
                                          update_gradients |
                                          update_normal_vectors |
                                          update_quadrature_points);

      Vector<double> local_rhs(dofs_per_cell);
      Vector<double> local_topography_adjoint(dofs_per_cell);
      Vector<double> local_cbf_mass_matrix(dofs_per_cell);
      Vector<double> local_cbf_state_adjoint(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<types::global_dof_index> face_dof_indices(this->get_fe().dofs_per_face);
      std::vector<double> topo_values(quadrature_formula_face.size());

      LinearAlgebra::BlockVector topography_adjoint(this->introspection().index_sets.system_partitioning,
                                                    this->introspection().index_sets.system_relevant_partitioning,
                                                    this->get_mpi_communicator());
      topography_adjoint = 0.0;

      LinearAlgebra::BlockVector cbf_mass_matrix(this->introspection().index_sets.system_partitioning,
                                                 this->introspection().index_sets.system_relevant_partitioning,
                                                 this->get_mpi_communicator());
      cbf_mass_matrix = 0.0;

      double local_top_area = 0.0;
      double local_surface_pressure_weight = 0.0;

      // Assemble J_h, the derivative of the objective with respect to
      // dynamic-topography values. For a misfit objective this seed is
      // h - h_target, not the computed topography h itself.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_face_values.reinit(cell, f);
                local_topography_adjoint = 0.0;

                MaterialModel::MaterialModelInputs<dim> in(fe_face_values,
                                                           cell,
                                                           this->introspection(),
                                                           this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> out(fe_face_values.n_quadrature_points,
                                                             this->n_compositional_fields());
                in.requested_properties = MaterialModel::MaterialProperties::density;
                this->get_material_model().evaluate(in, out);

                fe_face_values[this->introspection().extractors.temperature].get_function_values(dynamic_topography.topography_vector(),
                    topo_values);

                for (unsigned int q = 0; q < quadrature_formula_face.size(); ++q)
                  {
                    const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(fe_face_values.quadrature_point(q));
                    const double gravity_norm = gravity.norm();
                    const double density_contrast = out.densities[q] - dynamic_topography.get_density_above();

                    AssertThrow(std::abs(density_contrast) > std::numeric_limits<double>::min(),
                                ExcMessage("DynamicTopographyObjective encountered zero density contrast at the top boundary."));
                    AssertThrow(gravity_norm > std::numeric_limits<double>::min(),
                                ExcMessage("DynamicTopographyObjective encountered zero gravity magnitude at the top boundary."));

                    local_top_area += fe_face_values.JxW(q);
                    const double residual =
                      topo_values[q] - target_topography(fe_face_values.quadrature_point(q));
                    local_surface_pressure_weight +=
                      residual / (density_contrast * gravity_norm) * fe_face_values.JxW(q);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      local_topography_adjoint(i) +=
                        residual
                        * fe_face_values[this->introspection().extractors.temperature].value(i, q)
                        * fe_face_values.JxW(q);
                  }

                cell->distribute_local_to_global(local_topography_adjoint, topography_adjoint);

                fe_cbf_face_values.reinit(cell, f);
                local_cbf_mass_matrix = 0.0;
                for (unsigned int q = 0; q < cbf_face_quadrature_formula.size(); ++q)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      const Tensor<1, dim> phi_u =
                        fe_cbf_face_values[this->introspection().extractors.velocities].value(i, q);
                      local_cbf_mass_matrix(i) += phi_u * phi_u * fe_cbf_face_values.JxW(q);
                    }

                cell->distribute_local_to_global(local_cbf_mass_matrix, cbf_mass_matrix);
              }

      topography_adjoint.compress(VectorOperation::add);
      cbf_mass_matrix.compress(VectorOperation::add);

      const double top_area = Utilities::MPI::sum(local_top_area, this->get_mpi_communicator());
      const double surface_pressure_weight =
        Utilities::MPI::sum(local_surface_pressure_weight, this->get_mpi_communicator());

      AssertThrow(top_area > std::numeric_limits<double>::min(),
                  ExcMessage("DynamicTopographyObjective found zero top boundary area."));

      // Map J_h through the CBF postprocessing formula for dynamic
      // topography and then apply the transpose of the reaction operator
      // to obtain the adjoint RHS on velocity and pressure DoFs.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_values.reinit(cell);
                fe_support_values.reinit(cell, f);
                local_rhs = 0.0;
                local_cbf_state_adjoint = 0.0;

                MaterialModel::MaterialModelInputs<dim> support_in(fe_support_values,
                                                                   cell,
                                                                   this->introspection(),
                                                                   this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> support_out(fe_support_values.n_quadrature_points,
                                                                     this->n_compositional_fields());
                support_in.requested_properties = MaterialModel::MaterialProperties::density;
                this->get_material_model().evaluate(support_in, support_out);

                cell->get_dof_indices(local_dof_indices);
                cell->face(f)->get_dof_indices(face_dof_indices);

                for (unsigned int face_i = 0; face_i < face_dof_indices.size(); ++face_i)
                  {
                    const std::pair<unsigned int, unsigned int> component_index =
                      this->get_fe().face_system_to_component_index(face_i);
                    const unsigned int component = component_index.first;
                    const unsigned int support_index = component_index.second;

                    if (component == this->introspection().component_indices.temperature)
                      {
                        const Tensor<1,dim> normal = fe_support_values.normal_vector(support_index);
                        const Tensor<1,dim> gravity =
                          this->get_gravity_model().gravity_vector(fe_support_values.quadrature_point(support_index));
                        const double gravity_norm = gravity.norm();
                        const double density_contrast =
                          support_out.densities[support_index] - dynamic_topography.get_density_above();

                        AssertThrow(std::abs(density_contrast) > std::numeric_limits<double>::min(),
                                    ExcMessage("DynamicTopographyObjective encountered zero density contrast at the top boundary."));
                        AssertThrow(gravity_norm > std::numeric_limits<double>::min(),
                                    ExcMessage("DynamicTopographyObjective encountered zero gravity magnitude at the top boundary."));

                        const double coefficient =
                          topography_adjoint(face_dof_indices[face_i])
                          / (density_contrast * gravity_norm);

                        for (unsigned int velocity_face_i = 0; velocity_face_i < face_dof_indices.size(); ++velocity_face_i)
                          {
                            const std::pair<unsigned int, unsigned int> velocity_component_index =
                              this->get_fe().face_system_to_component_index(velocity_face_i);
                            const unsigned int velocity_component = velocity_component_index.first;

                            if (this->introspection().component_masks.velocities[velocity_component] == false)
                              continue;

                            const auto local_dof = std::find(local_dof_indices.begin(),
                                                             local_dof_indices.end(),
                                                             face_dof_indices[velocity_face_i]);
                            AssertThrow(local_dof != local_dof_indices.end(), ExcInternalError());

                            const double mass = cbf_mass_matrix(face_dof_indices[velocity_face_i]);
                            if (mass <= std::numeric_limits<double>::min())
                              continue;

                            const unsigned int local_index = local_dof - local_dof_indices.begin();
                            const Tensor<1,dim> phi_u =
                              fe_support_values[this->introspection().extractors.velocities].value(local_index,
                                                                                                   support_index);
                            local_cbf_state_adjoint(local_index) += coefficient * (phi_u * normal) / mass;
                          }
                      }
                  }

                // Apply the volume part of the transposed CBF reaction:
                // the velocity block contributes viscous stress terms and
                // the pressure block contributes the incompressibility term.
                MaterialModel::MaterialModelInputs<dim> volume_in(fe_values,
                                                                  cell,
                                                                  this->introspection(),
                                                                  this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> volume_out(fe_values.n_quadrature_points,
                                                                    this->n_compositional_fields());
                volume_in.requested_properties = MaterialModel::MaterialProperties::viscosity;
                this->get_material_model().evaluate(volume_in, volume_out);

                for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                  {
                    const double JxW = fe_values.JxW(q);

                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      if (local_cbf_state_adjoint(j) != 0.0)
                        {
                          const SymmetricTensor<2,dim> epsilon_phi_j =
                            fe_values[this->introspection().extractors.velocities].symmetric_gradient(j, q);
                          const double div_phi_j =
                            fe_values[this->introspection().extractors.velocities].divergence(j, q);

                          for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            {
                              const SymmetricTensor<2,dim> epsilon_phi_i =
                                fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);
                              const double pressure_phi_i =
                                fe_values[this->introspection().extractors.pressure].value(i, q);

                              local_rhs(i) += local_cbf_state_adjoint(j)
                                              * (2.0 * volume_out.viscosities[q] * (epsilon_phi_j * epsilon_phi_i)
                                                 - this->get_pressure_scaling() * div_phi_j * pressure_phi_i)
                                              * JxW;
                            }
                        }
                  }

                // Add the transpose contribution from subtracting the
                // surface-average pressure in the dynamic-topography formula.
                fe_face_values.reinit(cell, f);
                for (unsigned int q = 0; q < quadrature_formula_face.size(); ++q)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    local_rhs(i) += this->get_pressure_scaling()
                                    * surface_pressure_weight / top_area
                                    * fe_face_values[this->introspection().extractors.pressure].value(i, q)
                                    * fe_face_values.JxW(q);

                this->get_current_constraints().distribute_local_to_global(local_rhs,
                                                                           local_dof_indices,
                                                                           rhs);
              }

      rhs.compress(VectorOperation::add);

      if (this->get_parameters().pressure_normalization != "no")
        {
          // The forward pressure normalization removes a constant pressure
          // mode. Apply the transpose correction so the adjoint RHS is
          // consistent with that algebraic normalization.
          const unsigned int pressure_block_index =
            this->introspection().block_indices.pressure;
          const double inverse_pressure_scaling = 1.0 / this->get_pressure_scaling();
          const IndexSet local_pressure_dofs =
            rhs.block(pressure_block_index).locally_owned_elements();

          double local_pressure_constant_rhs_dot = 0.0;
          for (unsigned int k = 0; k < local_pressure_dofs.n_elements(); ++k)
            {
              const types::global_dof_index local_index =
                local_pressure_dofs.nth_index_in_set(k);
              local_pressure_constant_rhs_dot +=
                rhs.block(pressure_block_index)[local_index] * inverse_pressure_scaling;
            }

          const double pressure_constant_rhs_dot =
            Utilities::MPI::sum(local_pressure_constant_rhs_dot,
                                this->get_mpi_communicator());

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              for (const unsigned int f : cell->face_indices())
                if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                  {
                    fe_face_values.reinit(cell, f);
                    local_rhs = 0.0;
                    cell->get_dof_indices(local_dof_indices);

                    for (unsigned int q = 0; q < quadrature_formula_face.size(); ++q)
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        local_rhs(i) +=
                          -pressure_constant_rhs_dot
                          * this->get_pressure_scaling()
                          / top_area
                          * fe_face_values[this->introspection().extractors.pressure].value(i, q)
                          * fe_face_values.JxW(q);

                    this->get_current_constraints().distribute_local_to_global(local_rhs,
                                                                               local_dof_indices,
                                                                               rhs);
                  }

          rhs.compress(VectorOperation::add);
        }
    }



    template <int dim>
    void
    DynamicTopographyObjective<dim>::add_direct_kernel_contributions(const std::string &objective_name,
                                                                     KernelRepository<dim> &kernels) const
    {
      const Postprocess::DynamicTopography<dim> &dynamic_topography =
        this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::DynamicTopography<dim>>();

      AssertThrow(dynamic_topography.topography_vector().size() != 0,
                  ExcMessage("The dynamic topography postprocessor must run before calculating dynamic-topography adjoint surface kernels."));

      const unsigned int n_active_cells = this->get_triangulation().n_active_cells();
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);
      const QGauss<dim-1> face_quadrature_formula(quadrature_degree);
      const QGaussLobatto<dim-1> cbf_face_quadrature_formula(quadrature_degree);
      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       face_quadrature_formula,
                                       update_values |
                                       update_gradients |
                                       update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

      FEFaceValues<dim> fe_cbf_face_values(this->get_mapping(),
                                           this->get_fe(),
                                           cbf_face_quadrature_formula,
                                           update_values |
                                           update_JxW_values);

      const std::vector<Point<dim-1>> face_support_points =
        this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_face_support_points();
      const Quadrature<dim-1> support_quadrature(face_support_points);
      FEFaceValues<dim> fe_support_values(this->get_mapping(),
                                          this->get_fe(),
                                          support_quadrature,
                                          update_values |
                                          update_gradients |
                                          update_normal_vectors |
                                          update_quadrature_points);

      Vector<double> cell_volumes(n_active_cells);
      Vector<double> density_surface(n_active_cells);
      Vector<double> viscosity_surface(n_active_cells);
      Vector<double> local_topography_adjoint(dofs_per_cell);
      Vector<double> local_cbf_mass_matrix(dofs_per_cell);
      Vector<double> local_density_cbf_rhs(dofs_per_cell);
      Vector<double> local_viscosity_cbf_rhs(dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
      std::vector<types::global_dof_index> face_dof_indices(this->get_fe().dofs_per_face);
      std::vector<double> topo_values(face_quadrature_formula.size());
      std::vector<double> support_topo_values(support_quadrature.size());

      LinearAlgebra::BlockVector topography_adjoint(this->introspection().index_sets.system_partitioning,
                                                    this->introspection().index_sets.system_relevant_partitioning,
                                                    this->get_mpi_communicator());
      topography_adjoint = 0.0;

      LinearAlgebra::BlockVector cbf_mass_matrix(this->introspection().index_sets.system_partitioning,
                                                 this->introspection().index_sets.system_relevant_partitioning,
                                                 this->get_mpi_communicator());
      cbf_mass_matrix = 0.0;

      // Kernels are stored as cell averages, so first collect cell volumes
      // and divide accumulated cell-integrated contributions at the end.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            const unsigned int cell_index = cell->active_cell_index();
            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
              cell_volumes(cell_index) += fe_values.JxW(q);
          }

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_face_values.reinit(cell, f);
                local_topography_adjoint = 0.0;

                // Reassemble the direct-term seed J_h. It is the same
                // topography misfit used in evaluate() and in the adjoint RHS.
                fe_face_values[this->introspection().extractors.temperature].get_function_values(
                  dynamic_topography.topography_vector(),
                  topo_values);

                for (unsigned int q = 0; q < face_quadrature_formula.size(); ++q)
                  {
                    const double residual =
                      topo_values[q] - target_topography(fe_face_values.quadrature_point(q));
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      local_topography_adjoint(i) +=
                        residual
                        * fe_face_values[this->introspection().extractors.temperature].value(i, q)
                        * fe_face_values.JxW(q);
                  }

                cell->distribute_local_to_global(local_topography_adjoint, topography_adjoint);
              }

      topography_adjoint.compress(VectorOperation::add);

      // The CBF mass is diagonal at the support points and maps reaction
      // residuals to the normal-stress values used by dynamic topography.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_cbf_face_values.reinit(cell, f);
                local_cbf_mass_matrix = 0.0;

                for (unsigned int q = 0; q < cbf_face_quadrature_formula.size(); ++q)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    {
                      const Tensor<1, dim> phi_u =
                        fe_cbf_face_values[this->introspection().extractors.velocities].value(i, q);
                      local_cbf_mass_matrix(i) += phi_u * phi_u * fe_cbf_face_values.JxW(q);
                    }

                cell->distribute_local_to_global(local_cbf_mass_matrix, cbf_mass_matrix);
              }

      cbf_mass_matrix.compress(VectorOperation::add);

      std::map<types::global_dof_index, std::pair<unsigned int, double>> density_surface_by_topography_dof;
      // Direct density derivative of the denominator in
      // h = normal_stress / ((rho - rho_above) |g|). This term uses the
      // actual computed h, multiplied by the misfit seed J_h.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_support_values.reinit(cell, f);

                MaterialModel::MaterialModelInputs<dim> support_in(fe_support_values,
                                                                   cell,
                                                                   this->introspection(),
                                                                   this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> support_out(fe_support_values.n_quadrature_points,
                                                                     this->n_compositional_fields());
                support_in.requested_properties = MaterialModel::MaterialProperties::density;
                this->get_material_model().evaluate(support_in, support_out);

                fe_support_values[this->introspection().extractors.temperature].get_function_values(
                  dynamic_topography.topography_vector(),
                  support_topo_values);

                cell->face(f)->get_dof_indices(face_dof_indices);
                for (unsigned int face_i = 0; face_i < face_dof_indices.size(); ++face_i)
                  {
                    const std::pair<unsigned int, unsigned int> component_index =
                      this->get_fe().face_system_to_component_index(face_i);
                    const unsigned int component = component_index.first;
                    const unsigned int support_index = component_index.second;

                    if (component == this->introspection().component_indices.temperature)
                      {
                        const Tensor<1,dim> gravity =
                          this->get_gravity_model().gravity_vector(fe_support_values.quadrature_point(support_index));
                        const double gravity_norm = gravity.norm();
                        const double density_contrast =
                          support_out.densities[support_index] - dynamic_topography.get_density_above();

                        AssertThrow(std::abs(density_contrast) > std::numeric_limits<double>::min(),
                                    ExcMessage("DynamicTopographyObjective encountered zero density contrast at the top boundary."));
                        AssertThrow(gravity_norm > std::numeric_limits<double>::min(),
                                    ExcMessage("DynamicTopographyObjective encountered zero gravity magnitude at the top boundary."));

                        (void)gravity_norm;
                        density_surface_by_topography_dof[face_dof_indices[face_i]] =
                          std::make_pair(cell->active_cell_index(),
                                         -topography_adjoint(face_dof_indices[face_i])
                                         * support_topo_values[support_index]
                                         / density_contrast);
                      }
                  }
              }

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_values.reinit(cell);
                fe_support_values.reinit(cell, f);
                local_density_cbf_rhs = 0.0;
                local_viscosity_cbf_rhs = 0.0;

                // Direct material derivatives of the CBF reaction residual
                // with the Stokes state frozen: density changes buoyancy and
                // viscosity changes the viscous stress.
                MaterialModel::MaterialModelInputs<dim> volume_in(fe_values,
                                                                  cell,
                                                                  this->introspection(),
                                                                  this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> volume_out(fe_values.n_quadrature_points,
                                                                    this->n_compositional_fields());
                volume_in.requested_properties = MaterialModel::MaterialProperties::viscosity;
                this->get_material_model().evaluate(volume_in, volume_out);

                for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                  {
                    const Tensor<1, dim> gravity =
                      this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      {
                        const Tensor<1, dim> phi_u =
                          fe_values[this->introspection().extractors.velocities].value(i, q);
                        const SymmetricTensor<2, dim> epsilon_phi_u =
                          fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);

                        local_density_cbf_rhs(i) -=
                          (gravity * phi_u) * fe_values.JxW(q);
                        local_viscosity_cbf_rhs(i) +=
                          2.0 * (epsilon_phi_u * volume_in.strain_rate[q])
                          * fe_values.JxW(q);
                      }
                  }

                MaterialModel::MaterialModelInputs<dim> support_in(fe_support_values,
                                                                   cell,
                                                                   this->introspection(),
                                                                   this->get_solution());
                MaterialModel::MaterialModelOutputs<dim> support_out(fe_support_values.n_quadrature_points,
                                                                     this->n_compositional_fields());
                support_in.requested_properties = MaterialModel::MaterialProperties::density;
                this->get_material_model().evaluate(support_in, support_out);

                cell->get_dof_indices(local_dof_indices);
                cell->face(f)->get_dof_indices(face_dof_indices);

                for (unsigned int face_i = 0; face_i < face_dof_indices.size(); ++face_i)
                  {
                    const std::pair<unsigned int, unsigned int> component_index =
                      this->get_fe().face_system_to_component_index(face_i);
                    const unsigned int component = component_index.first;
                    const unsigned int support_index = component_index.second;

                    if (component == this->introspection().component_indices.temperature)
                      {
                        const Tensor<1, dim> normal = fe_support_values.normal_vector(support_index);
                        const Tensor<1, dim> gravity =
                          this->get_gravity_model().gravity_vector(fe_support_values.quadrature_point(support_index));
                        const double gravity_norm = gravity.norm();
                        const double density_contrast =
                          support_out.densities[support_index] - dynamic_topography.get_density_above();
                        double density_force_normal_stress_derivative = 0.0;
                        double normal_stress_derivative = 0.0;

                        AssertThrow(std::abs(density_contrast) > std::numeric_limits<double>::min(),
                                    ExcMessage("DynamicTopographyObjective encountered zero density contrast at the top boundary."));
                        AssertThrow(gravity_norm > std::numeric_limits<double>::min(),
                                    ExcMessage("DynamicTopographyObjective encountered zero gravity magnitude at the top boundary."));

                        for (unsigned int velocity_face_i = 0; velocity_face_i < face_dof_indices.size(); ++velocity_face_i)
                          {
                            const std::pair<unsigned int, unsigned int> velocity_component_index =
                              this->get_fe().face_system_to_component_index(velocity_face_i);
                            const unsigned int velocity_component = velocity_component_index.first;

                            if (this->introspection().component_masks.velocities[velocity_component] == false)
                              continue;

                            const auto local_dof = std::find(local_dof_indices.begin(),
                                                             local_dof_indices.end(),
                                                             face_dof_indices[velocity_face_i]);
                            AssertThrow(local_dof != local_dof_indices.end(), ExcInternalError());

                            const double mass = cbf_mass_matrix(face_dof_indices[velocity_face_i]);
                            if (mass <= std::numeric_limits<double>::min())
                              continue;

                            const unsigned int local_index = local_dof - local_dof_indices.begin();
                            const Tensor<1, dim> phi_u =
                              fe_support_values[this->introspection().extractors.velocities].value(local_index,
                                                                                                   support_index);
                            // Project the material derivative through the
                            // CBF mass inverse and the boundary normal.
                            density_force_normal_stress_derivative +=
                              (phi_u * normal)
                              * local_density_cbf_rhs(local_index)
                              / mass;
                            normal_stress_derivative +=
                              (phi_u * normal)
                              * local_viscosity_cbf_rhs(local_index)
                              / mass;
                          }

                        density_surface(cell->active_cell_index()) -=
                          topography_adjoint(face_dof_indices[face_i])
                          * density_force_normal_stress_derivative
                          / (density_contrast * gravity_norm);
                        viscosity_surface(cell->active_cell_index()) -=
                          topography_adjoint(face_dof_indices[face_i])
                          * normal_stress_derivative
                          / (density_contrast * gravity_norm);
                      }
                  }
              }

      // Add the denominator-only density contributions gathered by their
      // topography DoF before normalizing the cell kernels.
      for (const auto &entry : density_surface_by_topography_dof)
        density_surface(entry.second.first) += entry.second.second;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            const unsigned int cell_index = cell->active_cell_index();
            AssertThrow(cell_volumes(cell_index) > std::numeric_limits<double>::min(),
                        ExcMessage("DynamicTopographyObjective encountered a locally owned cell with zero volume."));

            density_surface(cell_index) /= cell_volumes(cell_index);
            viscosity_surface(cell_index) /= cell_volumes(cell_index);
          }

      kernels.add_contribution({objective_name,
                                "dynamic topography surface objective",
                                PhysicalProperty::density
                               },
                               density_surface);
      kernels.add_contribution({objective_name,
                                "dynamic topography surface objective",
                                PhysicalProperty::viscosity
                               },
                               viscosity_surface);
    }



    template <int dim>
    void
    DynamicTopographyObjective<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Adjoint");
      {
        prm.enter_subsection ("Dynamic topography");
        {
          prm.declare_entry ("Target topography model", "none",
                             Patterns::Selection ("none|function"),
                             "The model used to prescribe the target dynamic topography. "
                             "The dynamic topography objective evaluates one half of the "
                             "squared misfit between the computed and target topography. "
                             "A target topography model is required when this objective is used.");

          prm.enter_subsection ("Function");
          {
            prm.declare_entry ("Coordinate system", "spherical",
                               Patterns::Selection ("cartesian|spherical|depth"),
                               "A selection that determines the assumed coordinate "
                               "system for the target topography function variables. "
                               "Allowed values are `cartesian', `spherical', and `depth'. "
                               "`spherical' coordinates are interpreted as r,phi or "
                               "r,phi,theta in 2d/3d respectively with theta being the "
                               "polar angle.");

            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
            prm.declare_entry ("Variable names",
                               dim == 2 ? "r,phi" : "r,phi,theta",
                               Patterns::Anything(),
                               "Names of the coordinates used in the target topography function.");
            prm.declare_entry ("Function expression", "",
                               Patterns::Anything(),
                               "Target dynamic topography in meters.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DynamicTopographyObjective<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Adjoint");
      {
        prm.enter_subsection ("Dynamic topography");
        {
          target_topography_model = prm.get ("Target topography model");

          AssertThrow (target_topography_model != "none",
                       ExcMessage("Dynamic topography objective requires a target topography model. "
                                  "Set Adjoint/Dynamic topography/Target topography model."));

          if (target_topography_model == "function")
            {
              prm.enter_subsection ("Function");
              {
                target_coordinate_system =
                  Utilities::Coordinates::string_to_coordinate_system(prm.get ("Coordinate system"));
              }

              try
                {
                  target_topography_function.parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Adjoint.Dynamic topography.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get ("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }

              prm.leave_subsection();
            }
          else
            AssertThrow (false, ExcNotImplemented());
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

namespace aspect
{
  namespace Adjoint
  {
    ASPECT_REGISTER_ADJOINT_OBJECTIVE(DynamicTopographyObjective,
                                      "dynamic topography",
                                      "Objective functional for instantaneous dynamic-topography "
                                      "inversion. The objective evaluates one half of the "
                                      "squared surface dynamic topography misfit and assembles the "
                                      "corresponding top-boundary adjoint right hand side.")
  }
}
