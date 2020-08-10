/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/simulator.h>
#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/particle/world.h>

namespace aspect
{
  template <int dim>
  SimulatorAccess<dim>::SimulatorAccess ()
    :
    simulator (nullptr)
  {}


  template <int dim>
  SimulatorAccess<dim>::SimulatorAccess (const Simulator<dim> &simulator_object)
    :
    simulator (&simulator_object)
  {}


  template <int dim>
  SimulatorAccess<dim>::~SimulatorAccess ()
  {}



  template <int dim>
  void
  SimulatorAccess<dim>::initialize_simulator (const Simulator<dim> &simulator_object)
  {
    simulator = &simulator_object;
  }



  template <int dim>
  const Simulator<dim> &
  SimulatorAccess<dim>::get_simulator() const
  {
    return *simulator;
  }


  template <int dim>
  const Parameters<dim> &
  SimulatorAccess<dim>::get_parameters() const
  {
    return simulator->parameters;
  }


  template <int dim>
  SimulatorSignals<dim> &
  SimulatorAccess<dim>::get_signals() const
  {
    // Our reference to the Simulator is const, but we need to
    // be able to connect to the signals so a cast is required.
    return const_cast<SimulatorSignals<dim>&>(simulator->signals);
  }


  template <int dim>
  const Introspection<dim> &
  SimulatorAccess<dim>::introspection () const
  {
    return simulator->introspection;
  }


  template <int dim>
  MPI_Comm
  SimulatorAccess<dim>::get_mpi_communicator () const
  {
    return simulator->mpi_communicator;
  }

  template <int dim>
  TimerOutput &
  SimulatorAccess<dim>::get_computing_timer () const
  {
    return simulator->computing_timer;
  }

  template <int dim>
  const ConditionalOStream &
  SimulatorAccess<dim>::get_pcout () const
  {
    return simulator->pcout;
  }

  template <int dim>
  double SimulatorAccess<dim>::get_time () const
  {
    return simulator->time;
  }

  template <int dim>
  double SimulatorAccess<dim>::get_timestep () const
  {
    return simulator->time_step;
  }

  template <int dim>
  double SimulatorAccess<dim>::get_old_timestep () const
  {
    return simulator->old_time_step;
  }



  template <int dim>
  unsigned int SimulatorAccess<dim>::get_timestep_number () const
  {
    return simulator->timestep_number;
  }


  template <int dim>
  unsigned int SimulatorAccess<dim>::get_nonlinear_iteration () const
  {
    return simulator->nonlinear_iteration;
  }


  template <int dim>
  const parallel::distributed::Triangulation<dim> &
  SimulatorAccess<dim>::get_triangulation () const
  {
    return simulator->triangulation;
  }



  template <int dim>
  double
  SimulatorAccess<dim>::get_volume () const
  {
    return simulator->global_volume;
  }



  template <int dim>
  const Mapping<dim> &
  SimulatorAccess<dim>::get_mapping () const
  {
    return *(simulator->mapping);
  }



  template <int dim>
  std::string
  SimulatorAccess<dim>::get_output_directory () const
  {
    return simulator->parameters.output_directory;
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::convert_output_to_years () const
  {
    return simulator->parameters.convert_to_years;
  }


  template <int dim>
  unsigned int
  SimulatorAccess<dim>::get_pre_refinement_step () const
  {
    return simulator->pre_refinement_step;
  }


  template <int dim>
  unsigned int
  SimulatorAccess<dim>::n_compositional_fields () const
  {
    return simulator->introspection.n_compositional_fields;
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::include_adiabatic_heating () const
  {
    return simulator->heating_model_manager.adiabatic_heating_enabled();
  }

  template <int dim>
  bool
  SimulatorAccess<dim>::include_latent_heat () const
  {
    const std::vector<std::string> &heating_models = simulator->heating_model_manager.get_active_heating_model_names();
    return (std::find(heating_models.begin(), heating_models.end(), "latent heat") != heating_models.end());
  }

  template <int dim>
  bool
  SimulatorAccess<dim>::include_melt_transport () const
  {
    return simulator->parameters.include_melt_transport;
  }

  template <int dim>
  int
  SimulatorAccess<dim>::get_stokes_velocity_degree () const
  {
    return simulator->parameters.stokes_velocity_degree;
  }

  template <int dim>
  double
  SimulatorAccess<dim>::get_adiabatic_surface_temperature () const
  {
    return simulator->parameters.adiabatic_surface_temperature;
  }

  template <int dim>
  double
  SimulatorAccess<dim>::get_surface_pressure () const
  {
    return simulator->parameters.surface_pressure;
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_refinement_criteria (Vector<float> &estimated_error_per_cell) const
  {
    simulator->mesh_refinement_manager.execute (estimated_error_per_cell);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_artificial_viscosity (Vector<float> &viscosity_per_cell,
                                                  const bool skip_interior_cells) const
  {
    const typename Simulator<dim>::AdvectionField advection_field = Simulator<dim>::AdvectionField::temperature();
    simulator->get_artificial_viscosity(viscosity_per_cell, advection_field, skip_interior_cells);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_artificial_viscosity_composition (Vector<float> &viscosity_per_cell,
                                                              const unsigned int compositional_variable) const
  {
    const typename Simulator<dim>::AdvectionField advection_field = Simulator<dim>::AdvectionField::composition(compositional_variable);
    simulator->get_artificial_viscosity(viscosity_per_cell, advection_field);
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_current_linearization_point () const
  {
    return simulator->current_linearization_point;
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_solution () const
  {
    return simulator->solution;
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_old_solution () const
  {
    return simulator->old_solution;
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_old_old_solution () const
  {
    return simulator->old_old_solution;
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_reaction_vector () const
  {
    return simulator->operator_split_reaction_vector;
  }

  template <int dim>
  const LinearAlgebra::BlockVector &
  SimulatorAccess<dim>::get_mesh_velocity () const
  {
    Assert( simulator->parameters.mesh_deformation_enabled,
            ExcMessage("You cannot get the mesh velocity if mesh deformation is not enabled."));
    return simulator->mesh_deformation->mesh_velocity;
  }


  template <int dim>
  const DoFHandler<dim> &
  SimulatorAccess<dim>::get_dof_handler () const
  {
    return simulator->dof_handler;
  }



  template <int dim>
  const FiniteElement<dim> &
  SimulatorAccess<dim>::get_fe () const
  {
    Assert (simulator->dof_handler.n_dofs() > 0,
            ExcMessage("You are trying to access the FiniteElement before the DOFs have been "
                       "initialized. This may happen when accessing the Simulator from a plugin "
                       "that gets executed early in some cases (like material models) or from "
                       "an early point in the core code."));
    return simulator->dof_handler.get_fe();
  }

  template <int dim>
  const LinearAlgebra::BlockSparseMatrix &
  SimulatorAccess<dim>::get_system_matrix () const
  {
    return simulator->system_matrix;
  }

  template <int dim>
  const LinearAlgebra::BlockSparseMatrix &
  SimulatorAccess<dim>::get_system_preconditioner_matrix () const
  {
    return simulator->system_preconditioner_matrix;
  }

  template <int dim>
  const MaterialModel::Interface<dim> &
  SimulatorAccess<dim>::get_material_model () const
  {
    Assert (simulator->material_model.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->material_model.get();
  }


  template <int dim>
  void
  SimulatorAccess<dim>::compute_material_model_input_values (const LinearAlgebra::BlockVector                            &input_solution,
                                                             const FEValuesBase<dim,dim>                                 &input_finite_element_values,
                                                             const typename DoFHandler<dim>::active_cell_iterator        &cell,
                                                             const bool                                                   compute_strainrate,
                                                             MaterialModel::MaterialModelInputs<dim> &material_model_inputs) const
  {
    simulator->compute_material_model_input_values(input_solution,
                                                   input_finite_element_values,
                                                   cell,
                                                   compute_strainrate,
                                                   material_model_inputs);
  }



  template <int dim>
  const std::map<types::boundary_id,std::unique_ptr<BoundaryTraction::Interface<dim> > > &
  SimulatorAccess<dim>::get_boundary_traction () const
  {
    return simulator->boundary_traction;
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::has_boundary_temperature () const
  {
    return (get_boundary_temperature_manager().get_fixed_temperature_boundary_indicators().size() > 0);
  }



  template <int dim>
  const BoundaryTemperature::Interface<dim> &
  SimulatorAccess<dim>::get_boundary_temperature () const
  {
    Assert (get_boundary_temperature_manager().get_active_boundary_temperature_conditions().size() == 1,
            ExcMessage("You can only call this function if exactly one boundary temperature plugin is active."));
    return *(get_boundary_temperature_manager().get_active_boundary_temperature_conditions().front());
  }



  template <int dim>
  const BoundaryTemperature::Manager<dim> &
  SimulatorAccess<dim>::get_boundary_temperature_manager () const
  {
    return simulator->boundary_temperature_manager;
  }



  template <int dim>
  const BoundaryHeatFlux::Interface<dim> &
  SimulatorAccess<dim>::get_boundary_heat_flux () const
  {
    Assert (simulator->boundary_heat_flux.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->boundary_heat_flux.get();
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::has_boundary_composition () const
  {
    return (get_boundary_composition_manager().get_fixed_composition_boundary_indicators().size() > 0);
  }



  template <int dim>
  const BoundaryComposition::Interface<dim> &
  SimulatorAccess<dim>::get_boundary_composition () const
  {
    Assert (get_boundary_composition_manager().get_active_boundary_composition_conditions().size() == 1,
            ExcMessage("You can only call this function if exactly one boundary composition plugin is active."));
    return *(get_boundary_composition_manager().get_active_boundary_composition_conditions().front());
  }



  template <int dim>
  const BoundaryComposition::Manager<dim> &
  SimulatorAccess<dim>::get_boundary_composition_manager () const
  {
    return simulator->boundary_composition_manager;
  }



  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_temperature_boundary_indicators () const
  {
    return get_boundary_temperature_manager().get_fixed_temperature_boundary_indicators();
  }



  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_heat_flux_boundary_indicators () const
  {
    return simulator->parameters.fixed_heat_flux_boundary_indicators;
  }



  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_fixed_composition_boundary_indicators () const
  {
    return get_boundary_composition_manager().get_fixed_composition_boundary_indicators();
  }


  template <int dim>
  const std::set<types::boundary_id> &
  SimulatorAccess<dim>::get_mesh_deformation_boundary_indicators () const
  {
    Assert( simulator->parameters.mesh_deformation_enabled,
            ExcMessage("You cannot get the mesh deformation boundary indicators if mesh deformation is not enabled."));
    return simulator->mesh_deformation->get_active_mesh_deformation_boundary_indicators();
  }



  template <int dim>
  const BoundaryVelocity::Manager<dim> &
  SimulatorAccess<dim>::get_boundary_velocity_manager () const
  {
    return simulator->boundary_velocity_manager;
  }


  template <int dim>
  const InitialTopographyModel::Interface<dim> &
  SimulatorAccess<dim>::get_initial_topography_model () const
  {
    Assert (simulator->initial_topography_model.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->initial_topography_model.get();
  }


  template <int dim>
  const GeometryModel::Interface<dim> &
  SimulatorAccess<dim>::get_geometry_model () const
  {
    Assert (simulator->geometry_model.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->geometry_model.get();
  }

  template <int dim>
  const GravityModel::Interface<dim> &
  SimulatorAccess<dim>::get_gravity_model () const
  {
    Assert (simulator->gravity_model.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->gravity_model.get();
  }


  template <int dim>
  const AdiabaticConditions::Interface<dim> &
  SimulatorAccess<dim>::get_adiabatic_conditions () const
  {
    Assert (simulator->adiabatic_conditions.get() != nullptr,
            ExcMessage("You can not call this function if no such model is actually available."));
    return *simulator->adiabatic_conditions.get();
  }


  template <int dim>
  const InitialTemperature::Interface<dim> &
  SimulatorAccess<dim>::get_initial_temperature () const
  {
    Assert (get_initial_temperature_manager().get_active_initial_temperature_conditions().size() == 1,
            ExcMessage("You can only call this function if exactly one initial temperature plugin is active."));
    return *(get_initial_temperature_manager().get_active_initial_temperature_conditions().front());
  }


  template <int dim>
  const InitialTemperature::Manager<dim> &
  SimulatorAccess<dim>::get_initial_temperature_manager () const
  {
    return simulator->initial_temperature_manager;
  }


  template <int dim>
  const InitialComposition::Interface<dim> &
  SimulatorAccess<dim>::get_initial_composition () const
  {
    Assert (get_initial_composition_manager().get_active_initial_composition_conditions().size() == 1,
            ExcMessage("You can only call this function if only one initial composition plugin is active."));
    return *(get_initial_composition_manager().get_active_initial_composition_conditions().front());
  }


  template <int dim>
  const InitialComposition::Manager<dim> &
  SimulatorAccess<dim>::get_initial_composition_manager () const
  {
    return simulator->initial_composition_manager;
  }


  template <int dim>
  const HeatingModel::Manager<dim> &
  SimulatorAccess<dim>::get_heating_model_manager () const
  {
    return simulator->heating_model_manager;
  }

  template <int dim>
  const MeshRefinement::Manager<dim> &
  SimulatorAccess<dim>::get_mesh_refinement_manager () const
  {
    return simulator->mesh_refinement_manager;
  }

  template <int dim>
  const MeltHandler<dim> &
  SimulatorAccess<dim>::get_melt_handler () const
  {
    Assert (simulator->melt_handler.get() != nullptr,
            ExcMessage("You can not call this function if melt transport is not enabled."));
    return *(simulator->melt_handler);
  }

  template <int dim>
  const VolumeOfFluidHandler<dim> &
  SimulatorAccess<dim>::get_volume_of_fluid_handler () const
  {
    Assert (simulator->volume_of_fluid_handler.get() != nullptr,
            ExcMessage("You can not call this function if volume of fluid interface tracking is not enabled."));
    return *(simulator->volume_of_fluid_handler);
  }

  template <int dim>
  const NewtonHandler<dim> &
  SimulatorAccess<dim>::get_newton_handler () const
  {
    Assert (simulator->newton_handler.get() != nullptr,
            ExcMessage("You can not call this function if the Newton solver is not enabled."));
    return *(simulator->newton_handler);
  }


#ifdef ASPECT_WITH_WORLD_BUILDER
  template <int dim>
  const WorldBuilder::World &
  SimulatorAccess<dim>::get_world_builder () const
  {
    Assert (simulator->world_builder.get() != nullptr,
            ExcMessage("You can not call this function if the World Builder is not enabled. "
                       "Enable it by providing a path to a world builder file."));
    return *(simulator->world_builder);
  }
#endif


  template <int dim>
  const MeshDeformation::MeshDeformationHandler<dim> &
  SimulatorAccess<dim>::get_mesh_deformation_handler () const
  {
    Assert (simulator->mesh_deformation.get() != nullptr,
            ExcMessage("You cannot call this function if mesh deformation is not enabled."));

    return *(simulator->mesh_deformation);
  }

  template <int dim>
  void
  SimulatorAccess<dim>::get_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                                           const unsigned int                      q,
                                                           std::vector<double>                    &composition_values_at_q_point)
  {
    for (unsigned int k=0; k < composition_values_at_q_point.size(); ++k)
      composition_values_at_q_point[k] = composition_values[k][q];
  }


  template <int dim>
  TableHandler &
  SimulatorAccess<dim>::get_statistics_object () const
  {
    return const_cast<TableHandler &>(simulator->statistics);
  }

  template <int dim>
  const LateralAveraging<dim> &
  SimulatorAccess<dim>::get_lateral_averaging() const
  {
    return simulator->lateral_averaging;
  }

  template <int dim>
  const AffineConstraints<double> &
  SimulatorAccess<dim>::get_current_constraints() const
  {
    return simulator->current_constraints;
  }



  template <int dim>
  bool
  SimulatorAccess<dim>::simulator_is_past_initialization () const
  {
    return ((simulator != nullptr)
            &&
            (simulator->simulator_is_past_initialization == true));
  }


  template <int dim>
  double
  SimulatorAccess<dim>::get_pressure_scaling () const
  {
    return (simulator->pressure_scaling);
  }

  template <int dim>
  bool
  SimulatorAccess<dim>::pressure_rhs_needs_compatibility_modification () const
  {
    return simulator->do_pressure_rhs_compatibility_modification;
  }

  template <int dim>
  bool
  SimulatorAccess<dim>::model_has_prescribed_stokes_solution () const
  {
    return (simulator->prescribed_stokes_solution.get() != nullptr);
  }

  template <int dim>
  const Postprocess::Manager<dim> &
  SimulatorAccess<dim>::get_postprocess_manager() const
  {
    return simulator->postprocess_manager;
  }


  template <int dim>
  const Particle::World<dim> &
  SimulatorAccess<dim>::get_particle_world() const
  {
    Assert (simulator->particle_world.get() != nullptr,
            ExcMessage("You can not call this function if there is no particle world."));
    return *simulator->particle_world.get();
  }

  template <int dim>
  Particle::World<dim> &
  SimulatorAccess<dim>::get_particle_world()
  {
    Assert (simulator->particle_world.get() != nullptr,
            ExcMessage("You can not call this function if there is no particle world."));
    return *simulator->particle_world.get();
  }



  template <int dim>
  bool SimulatorAccess<dim>::is_stokes_matrix_free()
  {
    return (simulator->stokes_matrix_free ? true : false);
  }

  template <int dim>
  const StokesMatrixFreeHandler<dim> &
  SimulatorAccess<dim>::get_stokes_matrix_free () const
  {
    Assert (simulator->stokes_matrix_free.get() != nullptr,
            ExcMessage("You can not call this function if the matrix-free Stokes solver is not used."));
    return *(simulator->stokes_matrix_free);
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template class SimulatorAccess<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
