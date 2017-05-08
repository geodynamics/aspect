/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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
#include <aspect/assembly.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Assemblers
  {
    namespace
    {
      /* These functions implement a reduced form of the code from deal.II's TriaAccessor::measure().
       * In the 3d dG case, a call to face->measure() is not implemented for non-planar faces.
       * Since we only care about the scaling here, it is enough to have an approximation instead.
       * The 2d case remains unchanged.
       */
      double
      approximate_face_measure(const DoFHandler<2>::face_iterator &face)
      {
        return (face->vertex(0)-face->vertex(1)).norm();
      }

      double
      approximate_face_measure(const DoFHandler<3>::face_iterator &face)
      {
        const Tensor<1,3> v03 = face->vertex(3) - face->vertex(0);
        const Tensor<1,3> v12 = face->vertex(2) - face->vertex(1);
        const Tensor<1,3> twice_area = cross_product_3d(v03, v12);
        return 0.5 * twice_area.norm();
      }
    }

    template <int dim>
    void
    AdvectionAssembler<dim>::local_assemble_advection_system (const typename Simulator<dim>::AdvectionField &advection_field,
                                                              const double artificial_viscosity,
                                                              internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                                              internal::Assembly::CopyData::AdvectionSystem<dim> &data) const
    {
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      const bool   use_bdf2_scheme = (this->get_timestep_number() > 1);
      const double time_step = this->get_timestep();
      const double old_time_step = this->get_old_timestep();

      const double bdf2_factor = (use_bdf2_scheme)? ((2*time_step + old_time_step) /
                                                     (time_step + old_time_step)) : 1.0;

      const bool advection_field_is_temperature = advection_field.is_temperature();
      const unsigned int solution_component = advection_field.component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          // precompute the values of shape functions and their gradients.
          // We only need to look up values of shape functions if they
          // belong to 'our' component. They are zero otherwise anyway.
          // Note that we later only look at the values that we do set here.
          for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
            {
              if (fe.system_to_component_index(i).first == solution_component)
                {
                  scratch.grad_phi_field[i_advection] = scratch.finite_element_values[solution_field].gradient (i,q);
                  scratch.phi_field[i_advection]      = scratch.finite_element_values[solution_field].value (i,q);
                  ++i_advection;
                }
              ++i;
            }

          const double density_c_P              =
            ((advection_field_is_temperature)
             ?
             scratch.material_model_outputs.densities[q] *
             scratch.material_model_outputs.specific_heat[q]
             :
             1.0);

          Assert (density_c_P >= 0,
                  ExcMessage ("The product of density and c_P needs to be a "
                              "non-negative quantity."));

          const double conductivity =
            ((advection_field_is_temperature)
             ?
             scratch.material_model_outputs.thermal_conductivities[q]
             :
             0.0);
          const double latent_heat_LHS =
            ((advection_field_is_temperature)
             ?
             scratch.heating_model_outputs.lhs_latent_heat_terms[q]
             :
             0.0);
          Assert (density_c_P + latent_heat_LHS >= 0,
                  ExcMessage ("The sum of density times c_P and the latent heat contribution "
                              "to the left hand side needs to be a non-negative quantity."));

          const double gamma =
            ((advection_field_is_temperature)
             ?
             scratch.heating_model_outputs.heating_source_terms[q]
             :
             0.0);

          const double reaction_term =
            ((advection_field_is_temperature)
             ?
             0.0
             :
             scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]);

          const double field_term_for_rhs
            = (use_bdf2_scheme ?
               (scratch.old_field_values[q] *
                (1 + time_step/old_time_step)
                -
                scratch.old_old_field_values[q] *
                (time_step * time_step) /
                (old_time_step * (time_step + old_time_step)))
               :
               scratch.old_field_values[q])
              *
              (density_c_P + latent_heat_LHS);

          Tensor<1,dim> current_u = scratch.current_velocity_values[q];
          //Subtract off the mesh velocity for ALE corrections if necessary
          if (this->get_parameters().free_surface_enabled)
            current_u -= scratch.mesh_velocity_values[q];

          const double JxW = scratch.finite_element_values.JxW(q);

          // do the actual assembly. note that we only need to loop over the advection
          // shape functions because these are the only contributions we compute here
          for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
            {
              data.local_rhs(i)
              += (field_term_for_rhs * scratch.phi_field[i]
                  + time_step *
                  scratch.phi_field[i]
                  * gamma
                  + scratch.phi_field[i]
                  * reaction_term)
                 *
                 JxW;

              for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                {
                  data.local_matrix(i,j)
                  += (
                       (time_step * (conductivity + artificial_viscosity)
                        * (scratch.grad_phi_field[i] * scratch.grad_phi_field[j]))
                       + ((time_step * (scratch.phi_field[i] * (current_u * scratch.grad_phi_field[j])))
                          + (bdf2_factor * scratch.phi_field[i] * scratch.phi_field[j])) *
                       (density_c_P + latent_heat_LHS)
                     )
                     * JxW;
                }
            }
        }
    }



    template <int dim>
    std::vector<double>
    AdvectionAssembler<dim>::compute_advection_system_residual(const typename Simulator<dim>::AdvectionField     &advection_field,
                                                               internal::Assembly::Scratch::AdvectionSystem<dim> &scratch) const
    {
      const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
      std::vector<double> residuals(n_q_points);

      if (advection_field.is_temperature())
        this->get_heating_model_manager().evaluate(scratch.material_model_inputs,
                                                   scratch.material_model_outputs,
                                                   scratch.heating_model_outputs);

      for (unsigned int q=0; q < n_q_points; ++q)
        {
          const Tensor<1,dim> u = (scratch.old_velocity_values[q] +
                                   scratch.old_old_velocity_values[q]) / 2;

          const double dField_dt = (this->get_old_timestep() == 0.0) ? 0.0 :
                                   (
                                     ((scratch.old_field_values)[q] - (scratch.old_old_field_values)[q])
                                     / this->get_old_timestep());
          const double u_grad_field = u * (scratch.old_field_grads[q] +
                                           scratch.old_old_field_grads[q]) / 2;

          if (advection_field.is_temperature())
            {
              const double density       = scratch.material_model_outputs.densities[q];
              const double conductivity  = scratch.material_model_outputs.thermal_conductivities[q];
              const double c_P           = scratch.material_model_outputs.specific_heat[q];
              const double k_Delta_field = conductivity * (scratch.old_field_laplacians[q] +
                                                           scratch.old_old_field_laplacians[q]) / 2;

              const double gamma           = scratch.heating_model_outputs.heating_source_terms[q];
              const double latent_heat_LHS = scratch.heating_model_outputs.lhs_latent_heat_terms[q];

              residuals[q]
                = std::abs((density * c_P + latent_heat_LHS) * (dField_dt + u_grad_field) - k_Delta_field - gamma);
            }
          else
            {
              const double dreaction_term_dt = (this->get_old_timestep() == 0) ? 0.0 :
                                               scratch.material_model_outputs.reaction_terms[q][advection_field.compositional_variable]
                                               / this->get_old_timestep();

              residuals[q] = std::abs(dField_dt + u_grad_field - dreaction_term_dt);
            }
        }
      return residuals;
    }



    template <int dim>
    void
    AdvectionAssembler<dim>::local_assemble_discontinuous_advection_boundary_face_terms(const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int face_no,
        const typename Simulator<dim>::AdvectionField &advection_field,
        internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
        internal::Assembly::CopyData::AdvectionSystem<dim> &data) const
    {
      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int n_q_points    = scratch.face_finite_element_values->n_quadrature_points;

      const double time_step = this->get_timestep();

      if (!advection_field.is_discontinuous(introspection))
        return;

      // also have the number of dofs that correspond just to the element for
      // the system we are currently trying to assemble
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

      Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
      Assert (scratch.face_grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
      Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());

      const unsigned int solution_component = advection_field.component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      typename DoFHandler<dim>::face_iterator face = cell->face (face_no);

      if (((parameters.fixed_temperature_boundary_indicators.find(
              face->boundary_id()
            )
            != parameters.fixed_temperature_boundary_indicators.end())
           && (advection_field.is_temperature()))
          ||
          (( parameters.fixed_composition_boundary_indicators.find(
               face->boundary_id()
             )
             != parameters.fixed_composition_boundary_indicators.end())
           && (!advection_field.is_temperature())))
        {
          /*
           * We are in the case of a Dirichlet temperature or composition boundary.
           * In the temperature case, impose the Dirichlet value weakly using a matrix term
           * and RHS term. In the composition case, Dirichlet conditions can only be imposed
           * on inflow boundaries, and we only have the flow-dependent terms, so we only
           * assemble the corresponding flow-dependent and matrix and RHS terms
           * if we are on an inflow boundary.
           */

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              // precompute the values of shape functions and their gradients.
              // We only need to look up values of shape functions if they
              // belong to 'our' component. They are zero otherwise anyway.
              // Note that we later only look at the values that we do set here.
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      scratch.face_grad_phi_field[i_advection] = (*scratch.face_finite_element_values)[solution_field].gradient (i, q);
                      scratch.face_phi_field[i_advection]      = (*scratch.face_finite_element_values)[solution_field].value (i, q);
                      ++i_advection;
                    }
                  ++i;
                }

              const double density_c_P              =
                ((advection_field.is_temperature())
                 ?
                 scratch.face_material_model_outputs.densities[q] *
                 scratch.face_material_model_outputs.specific_heat[q]
                 :
                 1.0);

              Assert (density_c_P >= 0,
                      ExcMessage ("The product of density and c_P needs to be a "
                                  "non-negative quantity."));

              const double conductivity =
                ((advection_field.is_temperature())
                 ?
                 scratch.face_material_model_outputs.thermal_conductivities[q]
                 :
                 0.0);
              const double latent_heat_LHS =
                ((advection_field.is_temperature())
                 ?
                 scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                 :
                 0.0);
              Assert (density_c_P + latent_heat_LHS >= 0,
                      ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                  "to the left hand side needs to be a non-negative quantity."));

              const double penalty = (advection_field.is_temperature()
                                      ?
                                      parameters.discontinuous_penalty
                                      * parameters.temperature_degree
                                      * parameters.temperature_degree
                                      / approximate_face_measure(face)
                                      * conductivity
                                      / (density_c_P + latent_heat_LHS)
                                      :
                                      0.0);

              const double dirichlet_value = (advection_field.is_temperature()
                                              ?
                                              this->get_boundary_temperature().boundary_temperature(
                                                cell->face(face_no)->boundary_id(),
                                                scratch.face_finite_element_values->quadrature_point(q))
                                              :
                                              this->get_boundary_composition().boundary_composition(
                                                cell->face(face_no)->boundary_id(),
                                                scratch.face_finite_element_values->quadrature_point(q),
                                                advection_field.compositional_variable));

              Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
              //Subtract off the mesh velocity for ALE corrections if necessary
              if (parameters.free_surface_enabled)
                current_u -= scratch.face_mesh_velocity_values[q];

              /**
               * The discontinuous Galerkin method uses 2 types of jumps over edges:
               * undirected and directed jumps. Undirected jumps are dependent only
               * on the order of the numbering of cells. Directed jumps are dependent
               * on the direction of the flow. Thus the flow-dependent terms below are
               * only calculated if the edge is an inflow edge.
               */
              const bool inflow = ((current_u * scratch.face_finite_element_values->normal_vector(q)) < 0.);

              for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                {
                  data.local_rhs(i)
                  += (- time_step *  conductivity
                      * scratch.face_grad_phi_field[i]
                      * scratch.face_finite_element_values->normal_vector(q)
                      * dirichlet_value

                      + time_step
                      * (density_c_P + latent_heat_LHS)
                      * penalty
                      * scratch.face_phi_field[i]
                      * dirichlet_value

                      + (inflow
                         ?
                         - (density_c_P + latent_heat_LHS)
                         * time_step
                         * (current_u
                            * scratch.face_finite_element_values->normal_vector(q))
                         * dirichlet_value
                         * scratch.face_phi_field[i]
                         :
                         0.)
                     )
                     *
                     scratch.face_finite_element_values->JxW(q);

                  // local_matrix terms
                  for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                    {
                      data.local_matrix(i,j)
                      += (- time_step *  conductivity
                          * scratch.face_grad_phi_field[i]
                          * scratch.face_finite_element_values->normal_vector(q)
                          * scratch.face_phi_field[j]

                          - time_step *  conductivity
                          * scratch.face_grad_phi_field[j]
                          * scratch.face_finite_element_values->normal_vector(q)
                          * scratch.face_phi_field[i]

                          + time_step
                          * (density_c_P + latent_heat_LHS)
                          * penalty
                          * scratch.face_phi_field[i]
                          * scratch.face_phi_field[j]

                          + (inflow
                             ?
                             - (density_c_P + latent_heat_LHS)
                             * time_step
                             * (current_u
                                * scratch.face_finite_element_values->normal_vector(q))
                             * scratch.face_phi_field[i]
                             * scratch.face_phi_field[j]
                             :
                             0.)
                         )
                         * scratch.face_finite_element_values->JxW(q);
                    }
                }
            }
        }
      else if (cell->has_periodic_neighbor (face_no))
        {
          // Periodic temperature/composition term: consider the corresponding periodic faces as the case of interior faces
          this->local_assemble_discontinuous_advection_interior_face_terms(cell, face_no, advection_field, scratch, data);
        }
      else
        {
          //Neumann temperature term - no non-zero contribution as only homogeneous Neumann boundary conditions are implemented elsewhere for temperature
        }
    }



    template <int dim>
    void
    AdvectionAssembler<dim>::local_assemble_discontinuous_advection_interior_face_terms(const typename DoFHandler<dim>::active_cell_iterator &cell,
        const unsigned int face_no,
        const typename Simulator<dim>::AdvectionField &advection_field,
        internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
        internal::Assembly::CopyData::AdvectionSystem<dim> &data) const
    {
      const Parameters<dim> &parameters = this->get_parameters();
      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();

      const unsigned int n_q_points    = scratch.face_finite_element_values->n_quadrature_points;

      const double time_step = this->get_timestep();

      if (!advection_field.is_discontinuous(introspection))
        return;

      // also have the number of dofs that correspond just to the element for
      // the system we are currently trying to assemble
      const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int dofs_per_cell = fe.dofs_per_cell;

      Assert (advection_dofs_per_cell < scratch.face_finite_element_values->get_fe().dofs_per_cell, ExcInternalError());
      Assert (scratch.face_grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
      Assert (scratch.face_phi_field.size() == advection_dofs_per_cell, ExcInternalError());

      const unsigned int solution_component = advection_field.component_index(introspection);

      const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

      typename DoFHandler<dim>::face_iterator face = cell->face (face_no);

      //interior face or periodic face - no contribution on RHS

      const typename DoFHandler<dim>::cell_iterator
      neighbor = cell->neighbor_or_periodic_neighbor (face_no);
      //note: "neighbor" defined above is NOT active_cell_iterator, so this includes cells that are refined
      //for example: cell with periodic boundary.
      Assert (neighbor.state() == IteratorState::valid,
              ExcInternalError());
      const bool cell_has_periodic_neighbor = cell->has_periodic_neighbor (face_no);

      if (!(face->has_children()))
        {
          if (neighbor->level () == cell->level () &&
              neighbor->active() &&
              (((neighbor->is_locally_owned()) && (cell->index() < neighbor->index()))
               ||
               ((!neighbor->is_locally_owned()) && (cell->subdomain_id() < neighbor->subdomain_id()))))
            {
              Assert (cell->is_locally_owned(), ExcInternalError());
              //cell and neighbor are equal-sized, and cell has been chosen to assemble this face, so calculate from cell

              const unsigned int neighbor2 =
                (cell->has_periodic_neighbor(face_no)
                 ?
                 //how does the periodic neighbor talk about this cell?
                 cell->periodic_neighbor_of_periodic_neighbor( face_no )
                 :
                 //how does the neighbor talk about this cell?
                 cell->neighbor_of_neighbor(face_no));

              //set up neighbor values
              scratch.neighbor_face_finite_element_values->reinit (neighbor, neighbor2);

              this->compute_material_model_input_values (this->get_current_linearization_point(),
                                                         *scratch.neighbor_face_finite_element_values,
                                                         neighbor,
                                                         true,
                                                         scratch.neighbor_face_material_model_inputs);
              this->get_material_model().evaluate(scratch.neighbor_face_material_model_inputs,
                                                  scratch.neighbor_face_material_model_outputs);

              this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                         scratch.neighbor_face_material_model_outputs,
                                                         scratch.neighbor_face_heating_model_outputs);

              std::vector<types::global_dof_index> neighbor_dof_indices (dofs_per_cell);
              // get all dof indices on the neighbor, then extract those
              // that correspond to the solution_field we are interested in
              neighbor->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face][i_advection] = neighbor_dof_indices[i];
                      ++i_advection;
                    }
                  ++i;
                }
              data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face] = true;

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  // precompute the values of shape functions and their gradients.
                  // We only need to look up values of shape functions if they
                  // belong to 'our' component. They are zero otherwise anyway.
                  // Note that we later only look at the values that we do set here.
                  for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
                    {
                      if (fe.system_to_component_index(i).first == solution_component)
                        {
                          scratch.face_grad_phi_field[i_advection]          = (*scratch.face_finite_element_values)[solution_field].gradient (i, q);
                          scratch.face_phi_field[i_advection]               = (*scratch.face_finite_element_values)[solution_field].value (i, q);
                          scratch.neighbor_face_grad_phi_field[i_advection] = (*scratch.neighbor_face_finite_element_values)[solution_field].gradient (i, q);
                          scratch.neighbor_face_phi_field[i_advection]      = (*scratch.neighbor_face_finite_element_values)[solution_field].value (i, q);
                          ++i_advection;
                        }
                      ++i;
                    }

                  const double density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.densities[q] *
                     scratch.face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  Assert (density_c_P >= 0,
                          ExcMessage ("The product of density and c_P needs to be a "
                                      "non-negative quantity."));

                  const double conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  Assert (density_c_P + latent_heat_LHS >= 0,
                          ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                      "to the left hand side needs to be a non-negative quantity."));

                  const double penalty = (advection_field.is_temperature()
                                          ?
                                          parameters.discontinuous_penalty
                                          * parameters.temperature_degree
                                          * parameters.temperature_degree
                                          / approximate_face_measure(face)
                                          * conductivity
                                          / (density_c_P + latent_heat_LHS)
                                          :
                                          0.0);

                  Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                  //Subtract off the mesh velocity for ALE corrections if necessary
                  if (parameters.free_surface_enabled)
                    current_u -= scratch.face_mesh_velocity_values[q];

                  const double neighbor_density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.densities[q] *
                     scratch.neighbor_face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  Assert (neighbor_density_c_P >= 0,
                          ExcMessage ("The product of density and c_P on the neighbor needs to be a "
                                      "non-negative quantity."));

                  const double neighbor_conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double neighbor_latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  Assert (neighbor_density_c_P + neighbor_latent_heat_LHS >= 0,
                          ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                      "to the left hand side on the neighbor needs to be a non-negative quantity."));

                  const double neighbor_penalty = (advection_field.is_temperature()
                                                   ?
                                                   parameters.discontinuous_penalty
                                                   * parameters.temperature_degree
                                                   * parameters.temperature_degree
                                                   / approximate_face_measure(neighbor->face(neighbor2))
                                                   * neighbor_conductivity
                                                   / (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                                   :
                                                   0.0);

                  const double max_penalty = std::max(penalty, neighbor_penalty);

                  const double max_density_c_P_and_latent_heat =
                    std::max(density_c_P + latent_heat_LHS,
                             neighbor_density_c_P + neighbor_latent_heat_LHS);

                  Assert (numbers::is_finite(max_density_c_P_and_latent_heat),
                          ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a finite quantity."));
                  Assert (max_density_c_P_and_latent_heat >= 0,
                          ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a "
                                      "non-negative quantity."));

                  /**
                   * The discontinuous Galerkin method uses 2 types of jumps over edges:
                   * undirected and directed jumps. Undirected jumps are dependent only
                   * on the order of the numbering of cells. Directed jumps are dependent
                   * on the direction of the flow. Thus the flow-dependent terms below are
                   * only calculated if the edge is an inflow edge.
                   */
                  const bool inflow = ((current_u * scratch.face_finite_element_values->normal_vector(q)) < 0.);

                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    {
                      for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j)
                          += (- 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[i]
                              * scratch.face_phi_field[j]

                              - (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);

                          data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face](i,j)
                          += (- 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[j]
                              * scratch.face_phi_field[i]

                              + (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);

                          data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face](i,j)
                          += (+ 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              - 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[j]
                              * scratch.neighbor_face_phi_field[i]

                              - (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);

                          data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face](i,j)
                          += (+ 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              + 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.face_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[i]
                              * scratch.neighbor_face_phi_field[j]

                              + (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.face_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.face_finite_element_values->JxW(q);
                        }
                    }
                }
            }
          else
            {
              /* neighbor is taking responsibility for assembly of this face, because
               * either (1) neighbor is coarser, or
               *        (2) neighbor is equally-sized and
               *           (a) neighbor is on a different subdomain, with lower subdmain_id(), or
               *           (b) neighbor is on the same subdomain and has lower index().
              */
            }
        }
      else //face->has_children(), so always assemble from here.
        {
          const unsigned int neighbor2 =
            (cell_has_periodic_neighbor
             ?
             cell->periodic_neighbor_face_no(face_no)
             :
             cell->neighbor_face_no(face_no));

          //loop over subfaces
          for (unsigned int subface_no=0; subface_no<face->number_of_children(); ++subface_no)
            {
              const typename DoFHandler<dim>::active_cell_iterator neighbor_child
                = ( cell_has_periodic_neighbor
                    ?
                    cell->periodic_neighbor_child_on_subface(face_no,subface_no)
                    :
                    cell->neighbor_child_on_subface (face_no, subface_no));

              //set up subface values
              scratch.subface_finite_element_values->reinit (cell, face_no, subface_no);

              //subface->face
              (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_current_linearization_point(),
                  scratch.face_current_velocity_values);

              //get the mesh velocity, as we need to subtract it off of the advection systems
              if (parameters.free_surface_enabled)
                (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_mesh_velocity(),
                    scratch.face_mesh_velocity_values);

              //get the mesh velocity, as we need to subtract it off of the advection systems
              if (parameters.free_surface_enabled)
                (*scratch.subface_finite_element_values)[introspection.extractors.velocities].get_function_values(this->get_mesh_velocity(),
                    scratch.face_mesh_velocity_values);

              this->compute_material_model_input_values (this->get_current_linearization_point(),
                                                         *scratch.subface_finite_element_values,
                                                         cell,
                                                         true,
                                                         scratch.face_material_model_inputs);
              this->get_material_model().evaluate(scratch.face_material_model_inputs,
                                                  scratch.face_material_model_outputs);

              this->get_heating_model_manager().evaluate(scratch.face_material_model_inputs,
                                                         scratch.face_material_model_outputs,
                                                         scratch.face_heating_model_outputs);

              //set up neighbor values
              scratch.neighbor_face_finite_element_values->reinit (neighbor_child, neighbor2);

              this->compute_material_model_input_values (this->get_current_linearization_point(),
                                                         *scratch.neighbor_face_finite_element_values,
                                                         neighbor_child,
                                                         true,
                                                         scratch.neighbor_face_material_model_inputs);
              this->get_material_model().evaluate(scratch.neighbor_face_material_model_inputs,
                                                  scratch.neighbor_face_material_model_outputs);

              this->get_heating_model_manager().evaluate(scratch.neighbor_face_material_model_inputs,
                                                         scratch.neighbor_face_material_model_outputs,
                                                         scratch.neighbor_face_heating_model_outputs);

              std::vector<types::global_dof_index> neighbor_dof_indices (fe.dofs_per_cell);
              // get all dof indices on the neighbor, then extract those
              // that correspond to the solution_field we are interested in
              neighbor_child->get_dof_indices (neighbor_dof_indices);
              for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell;/*increment at end of loop*/)
                {
                  if (fe.system_to_component_index(i).first == solution_component)
                    {
                      data.neighbor_dof_indices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no][i_advection] = neighbor_dof_indices[i];
                      ++i_advection;
                    }
                  ++i;
                }
              data.assembled_matrices[face_no * GeometryInfo<dim>::max_children_per_face + subface_no] = true;

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  // precompute the values of shape functions and their gradients.
                  // We only need to look up values of shape functions if they
                  // belong to 'our' component. They are zero otherwise anyway.
                  // Note that we later only look at the values that we do set here.
                  for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
                    {
                      if (fe.system_to_component_index(i).first == solution_component)
                        {
                          scratch.face_grad_phi_field[i_advection]          = (*scratch.subface_finite_element_values)[solution_field].gradient (i, q);
                          scratch.face_phi_field[i_advection]               = (*scratch.subface_finite_element_values)[solution_field].value (i, q);
                          scratch.neighbor_face_grad_phi_field[i_advection] = (*scratch.neighbor_face_finite_element_values)[solution_field].gradient (i, q);
                          scratch.neighbor_face_phi_field[i_advection]      = (*scratch.neighbor_face_finite_element_values)[solution_field].value (i, q);
                          ++i_advection;
                        }
                      ++i;
                    }

                  const double density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.densities[q] *
                     scratch.face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  Assert (density_c_P >= 0,
                          ExcMessage ("The product of density and c_P needs to be a "
                                      "non-negative quantity."));

                  const double conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  Assert (density_c_P + latent_heat_LHS >= 0,
                          ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                      "to the left hand side needs to be a non-negative quantity."));

                  const double penalty = (advection_field.is_temperature()
                                          ?
                                          parameters.discontinuous_penalty
                                          * parameters.temperature_degree
                                          * parameters.temperature_degree
                                          / approximate_face_measure(face)
                                          * conductivity
                                          / (density_c_P + latent_heat_LHS)
                                          :
                                          0.0);

                  Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];
                  //Subtract off the mesh velocity for ALE corrections if necessary
                  if (parameters.free_surface_enabled)
                    current_u -= scratch.face_mesh_velocity_values[q];

                  const double neighbor_density_c_P              =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.densities[q] *
                     scratch.neighbor_face_material_model_outputs.specific_heat[q]
                     :
                     1.0);

                  Assert (neighbor_density_c_P >= 0,
                          ExcMessage ("The product of density and c_P on the neighbor needs to be a "
                                      "non-negative quantity."));

                  const double neighbor_conductivity =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_material_model_outputs.thermal_conductivities[q]
                     :
                     0.0);
                  const double neighbor_latent_heat_LHS =
                    ((advection_field.is_temperature())
                     ?
                     scratch.neighbor_face_heating_model_outputs.lhs_latent_heat_terms[q]
                     :
                     0.0);
                  Assert (neighbor_density_c_P + neighbor_latent_heat_LHS >= 0,
                          ExcMessage ("The sum of density times c_P and the latent heat contribution "
                                      "to the left hand side on the neighbor needs to be a non-negative quantity."));

                  const double neighbor_penalty = (advection_field.is_temperature()
                                                   ?
                                                   parameters.discontinuous_penalty
                                                   * parameters.temperature_degree
                                                   * parameters.temperature_degree
                                                   / approximate_face_measure(neighbor_child->face(neighbor2))
                                                   * neighbor_conductivity
                                                   / (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                                   :
                                                   0.0);

                  const double max_penalty = std::max(penalty, neighbor_penalty);

                  const double max_density_c_P_and_latent_heat =
                    std::max(density_c_P + latent_heat_LHS,
                             neighbor_density_c_P + neighbor_latent_heat_LHS);

                  Assert (numbers::is_finite(max_density_c_P_and_latent_heat),
                          ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a finite quantity."));
                  Assert (max_density_c_P_and_latent_heat >= 0,
                          ExcMessage ("The maximum product of density and c_P plus latent heat LHS on the neighbor needs to be a "
                                      "non-negative quantity."));

                  /**
                   * The discontinuous Galerkin method uses 2 types of jumps over edges:
                   * undirected and directed jumps. Undirected jumps are dependent only
                   * on the order of the numbering of cells. Directed jumps are dependent
                   * on the direction of the flow. Thus the flow-dependent terms below are
                   * only calculated if the edge is an inflow edge.
                   */
                  const bool inflow = ((current_u * scratch.subface_finite_element_values->normal_vector(q)) < 0.);

                  for (unsigned int i=0; i<advection_dofs_per_cell; ++i)
                    {
                      for (unsigned int j=0; j<advection_dofs_per_cell; ++j)
                        {
                          data.local_matrix(i,j)
                          += (- 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[i]
                              * scratch.face_phi_field[j]

                              - (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);

                          data.local_matrices_int_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                          += (- 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[i]

                              + 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[j]
                              * scratch.face_phi_field[i]

                              + (inflow
                                 ?
                                 (density_c_P + latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);

                          data.local_matrices_ext_int[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                          += (+ 0.5 * time_step * conductivity
                              * scratch.face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              - 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.face_phi_field[j]

                              - time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.face_phi_field[j]
                              * scratch.neighbor_face_phi_field[i]

                              - (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);

                          data.local_matrices_ext_ext[face_no * GeometryInfo<dim>::max_children_per_face + subface_no](i,j)
                          += (+ 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[i]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[j]

                              + 0.5 * time_step * neighbor_conductivity
                              * scratch.neighbor_face_grad_phi_field[j]
                              * scratch.subface_finite_element_values->normal_vector(q)
                              * scratch.neighbor_face_phi_field[i]

                              + time_step
                              * max_density_c_P_and_latent_heat
                              * max_penalty
                              * scratch.neighbor_face_phi_field[i]
                              * scratch.neighbor_face_phi_field[j]

                              + (!inflow
                                 ?
                                 (neighbor_density_c_P + neighbor_latent_heat_LHS)
                                 * time_step
                                 * (current_u
                                    * scratch.subface_finite_element_values->normal_vector(q))
                                 * scratch.neighbor_face_phi_field[i]
                                 * scratch.neighbor_face_phi_field[j]
                                 :
                                 0.)
                             )
                             * scratch.subface_finite_element_values->JxW(q);
                        }
                    }
                }
            }
        }
    }
  }
} // namespace aspect

// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class \
  AdvectionAssembler<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
