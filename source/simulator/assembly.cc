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
#include <aspect/utilities.h>
#include <aspect/compat.h>
#include <aspect/simulator_access.h>
#include <aspect/citation_info.h>

#include <aspect/simulator/assemblers/interface.h>
#include <aspect/melt.h>
#include <aspect/newton.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator/assemblers/advection.h>

#include <aspect/stokes_matrix_free.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>

#include <limits>


namespace aspect
{
  template <int dim>
  void
  Simulator<dim>::
  compute_material_model_input_values (const LinearAlgebra::BlockVector                            &input_solution,
                                       const FEValuesBase<dim>                                     &input_finite_element_values,
                                       const typename DoFHandler<dim>::active_cell_iterator        &cell,
                                       const bool                                                   compute_strainrate,
                                       MaterialModel::MaterialModelInputs<dim>                     &material_model_inputs) const
  {
    const unsigned int n_q_points = material_model_inputs.temperature.size();

    material_model_inputs.position = input_finite_element_values.get_quadrature_points();

    input_finite_element_values[introspection.extractors.temperature].get_function_values (input_solution,
        material_model_inputs.temperature);
    input_finite_element_values[introspection.extractors.pressure].get_function_values(input_solution,
        material_model_inputs.pressure);
    input_finite_element_values[introspection.extractors.velocities].get_function_values(input_solution,
        material_model_inputs.velocity);
    input_finite_element_values[introspection.extractors.pressure].get_function_gradients (input_solution,
        material_model_inputs.pressure_gradient);

    // only the viscosity in the material can depend on the strain_rate
    // if this is not needed, we can save some time here. By setting the
    // length of the strain_rate vector to 0, we signal to evaluate()
    // that we do not need to access the viscosity.
    if (compute_strainrate)
      input_finite_element_values[introspection.extractors.velocities].get_function_symmetric_gradients(input_solution,
          material_model_inputs.strain_rate);
    else
      material_model_inputs.strain_rate.resize(0);

    // the values of the compositional fields are stored as block vectors for each field
    // we have to extract them in this structure
    std::vector<std::vector<double> > composition_values (introspection.n_compositional_fields,
                                                          std::vector<double> (n_q_points));

    for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
      input_finite_element_values[introspection.extractors.compositional_fields[c]].get_function_values(input_solution,
          composition_values[c]);

    // then we copy these values to exchange the inner and outer vector, because for the material
    // model we need a vector with values of all the compositional fields for every quadrature point
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<introspection.n_compositional_fields; ++c)
        material_model_inputs.composition[q][c] = composition_values[c][q];

    material_model_inputs.current_cell = cell;
  }



  namespace
  {
    // This function initializes the simulator access for all assemblers
    // inside of the assemblers parameter. This is just a shortcut to save some
    // lines in set_assemblers(), where this operation appears multiple times.
    template <int dim, class AssemblerType>
    void
    initialize_simulator(const Simulator<dim> &simulator,
                         std::vector<std::unique_ptr<AssemblerType> > &assemblers)
    {
      for (unsigned int i=0; i<assemblers.size(); ++i)
        if (SimulatorAccess<dim> *p = dynamic_cast<SimulatorAccess<dim>* >(assemblers[i].get()))
          p->initialize_simulator(simulator);
    }
  }

  template <int dim>
  void
  Simulator<dim>::
  set_stokes_assemblers()
  {
    assemblers->stokes_preconditioner.push_back(std_cxx14::make_unique<aspect::Assemblers::StokesPreconditioner<dim> >());
    assemblers->stokes_system.push_back(std_cxx14::make_unique<aspect::Assemblers::StokesIncompressibleTerms<dim> >());

    if (material_model->is_compressible())
      {
        // The compressible part of the preconditioner is only necessary if we use the simplified A block
        if (parameters.use_full_A_block_preconditioner == false)
          assemblers->stokes_preconditioner.push_back(
            std_cxx14::make_unique<aspect::Assemblers::StokesCompressiblePreconditioner<dim> >());

        assemblers->stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesCompressibleStrainRateViscosityTerm<dim> >());
      }

    if (parameters.formulation_mass_conservation ==
        Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile)
      {
        assemblers->stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesImplicitReferenceDensityCompressibilityTerm<dim> >());
      }
    else if (parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::reference_density_profile)
      {
        assemblers->stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesReferenceDensityCompressibilityTerm<dim> >());
      }
    else if (parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::incompressible)
      {
        // do nothing, because we assembled div u =0 above already
      }
    else if (parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::isentropic_compression)
      {
        assemblers->stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesIsentropicCompressionTerm<dim> >());
      }
    else if (parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::hydrostatic_compression)
      {
        assemblers->stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesHydrostaticCompressionTerm<dim> >());
      }
    else if (parameters.formulation_mass_conservation ==
             Parameters<dim>::Formulation::MassConservation::projected_density_field)
      {
        CitationInfo::add("pda");
        assemblers->stokes_system.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesProjectedDensityFieldTerm<dim> >());
      }
    else
      AssertThrow(false,
                  ExcMessage("Unknown mass conservation equation approximation. There is no assembler"
                             " defined that handles this formulation."));

    // add the terms for traction boundary conditions
    if (!boundary_traction.empty())
      {
        assemblers->stokes_system_on_boundary_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::StokesBoundaryTraction<dim> >());
      }

    // add the terms necessary to normalize the pressure
    if (do_pressure_rhs_compatibility_modification)
      assemblers->stokes_system.push_back(
        std_cxx14::make_unique<aspect::Assemblers::StokesPressureRHSCompatibilityModification<dim> >());

  }

  template <int dim>
  void
  Simulator<dim>::
  set_advection_assemblers()
  {
    assemblers->advection_system.push_back(
      std_cxx14::make_unique<aspect::Assemblers::AdvectionSystem<dim> >());

    // add the diffusion assemblers if we have fields that use this method
    if (std::find(parameters.compositional_field_methods.begin(), parameters.compositional_field_methods.end(),
                  Parameters<dim>::AdvectionFieldMethod::prescribed_field_with_diffusion) != parameters.compositional_field_methods.end())
      assemblers->advection_system.push_back(
        std_cxx14::make_unique<aspect::Assemblers::DiffusionSystem<dim> >());

    if (parameters.use_discontinuous_temperature_discretization ||
        parameters.use_discontinuous_composition_discretization)
      {
        // TODO: This currently does not work in parallel, because the sparsity
        // pattern of the matrix does not seem to know about flux terms
        // across periodic faces of different levels. Fix this.
        AssertThrow(geometry_model->get_periodic_boundary_pairs().size() == 0 ||
                    Utilities::MPI::n_mpi_processes(mpi_communicator) == 1 ||
                    (parameters.initial_adaptive_refinement == 0 &&
                     parameters.adaptive_refinement_interval == 0),
                    ExcMessage("Combining discontinuous elements with periodic boundaries and "
                               "adaptive mesh refinement in parallel models is currently not supported. "
                               "Please switch off any of those options or run on a single process."));

        assemblers->advection_system_on_boundary_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::AdvectionSystemBoundaryFace<dim> >());

        assemblers->advection_system_on_interior_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::AdvectionSystemInteriorFace<dim> >());
      }

    if (parameters.fixed_heat_flux_boundary_indicators.size() != 0)
      {
        assemblers->advection_system_on_boundary_face.push_back(
          std_cxx14::make_unique<aspect::Assemblers::AdvectionSystemBoundaryHeatFlux<dim> >());
      }

    if (parameters.use_discontinuous_temperature_discretization
        || parameters.fixed_heat_flux_boundary_indicators.size() != 0)
      {
        assemblers->advection_system_assembler_on_face_properties[0].need_face_material_model_data = true;
        assemblers->advection_system_assembler_on_face_properties[0].need_face_finite_element_evaluation = true;
      }

    if (parameters.use_discontinuous_composition_discretization)
      {
        for (unsigned int i = 1; i<=introspection.n_compositional_fields; ++i)
          {
            assemblers->advection_system_assembler_on_face_properties[i].need_face_material_model_data = true;
            assemblers->advection_system_assembler_on_face_properties[i].need_face_finite_element_evaluation = true;
          }
      }
  }

  template <int dim>
  void
  Simulator<dim>::
  set_assemblers ()
  {
    // first let the manager delete all existing assemblers:
    assemblers->reset();

    assemblers->advection_system_assembler_properties.resize(1+introspection.n_compositional_fields);
    assemblers->advection_system_assembler_on_face_properties.resize(1+introspection.n_compositional_fields);

    // Set up the set of assemblers depending on the mode of equations
    if (!parameters.include_melt_transport
        && !assemble_newton_stokes_system)
      {
        set_advection_assemblers();
        set_stokes_assemblers();
      }
    else if (parameters.include_melt_transport
             && !assemble_newton_stokes_system)
      {
        melt_handler->set_assemblers(*assemblers);
      }
    else if (!parameters.include_melt_transport
             && assemble_newton_stokes_system)
      {
        set_advection_assemblers();
        newton_handler->set_assemblers(*assemblers);
      }
    else if (parameters.include_melt_transport
             && assemble_newton_stokes_system)
      {
        AssertThrow (false,
                     ExcMessage("The melt implementation does not support the Newton solver at the moment."));
      }
    else
      AssertThrow(false,ExcInternalError());

    // allow other assemblers to add themselves or modify the existing ones by firing the signal
    this->signals.set_assemblers(*this, *assemblers);

    // ensure that all assembler objects have access to the SimulatorAccess
    // base class
    initialize_simulator(*this,assemblers->stokes_preconditioner);
    initialize_simulator(*this,assemblers->stokes_system);
    initialize_simulator(*this,assemblers->stokes_system_on_boundary_face);
    initialize_simulator(*this,assemblers->advection_system);
    initialize_simulator(*this,assemblers->advection_system_on_boundary_face);
    initialize_simulator(*this,assemblers->advection_system_on_interior_face);
  }


  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                        internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    // First get all dof indices of the current cell, then extract those
    // that correspond to the Stokes system we are interested in.
    // Note that assemblers below can modify this list of dofs, if they in fact
    // assemble a different system than the standard Stokes system (e.g. in
    // models with melt transport).

    cell->get_dof_indices (scratch.local_dof_indices);
    data.extract_stokes_dof_indices(scratch.local_dof_indices, introspection, finite_element);

    // Prepare the data structures for assembly
    scratch.reinit(cell);
    data.local_matrix = 0;

    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         cell,
                                         true,
                                         scratch.material_model_inputs);

    for (unsigned int i=0; i<assemblers->stokes_preconditioner.size(); ++i)
      assemblers->stokes_preconditioner[i]->create_additional_material_model_outputs(scratch.material_model_outputs);

    material_model->evaluate(scratch.material_model_inputs,
                             scratch.material_model_outputs);
    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.finite_element_values.get_quadrature(),
                                               scratch.finite_element_values.get_mapping(),
                                               scratch.material_model_outputs);

    for (unsigned int i=0; i<assemblers->stokes_preconditioner.size(); ++i)
      assemblers->stokes_preconditioner[i]->execute(scratch,data);
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_preconditioner (const internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_dof_indices,
                                                    system_preconditioner_matrix);
  }



  template <int dim>
  void
  Simulator<dim>::assemble_stokes_preconditioner ()
  {
    if (stokes_matrix_free)
      return;

    system_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    // determine which update flags to use for the cell integrals
    const UpdateFlags cell_update_flags
      = ((update_JxW_values |
          update_values |
          update_gradients |
          update_quadrature_points)
         |
         assemblers->stokes_preconditioner_assembler_properties.needed_update_flags);

    unsigned int stokes_dofs_per_cell = dim * finite_element.base_element(introspection.base_elements.velocities).dofs_per_cell
                                        + finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell;

    if (parameters.include_melt_transport)
      stokes_dofs_per_cell += finite_element.base_element(introspection.variable("compaction pressure").base_index).dofs_per_cell;

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                      internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
    {
      this->local_assemble_stokes_preconditioner(cell, scratch, data);
    };

    auto copier = [&](const internal::Assembly::CopyData::StokesPreconditioner<dim> &data)
    {
      this->copy_local_to_global_stokes_preconditioner(data);
    };

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         worker,
         copier,
         internal::Assembly::Scratch::
         StokesPreconditioner<dim> (finite_element, quadrature_formula,
                                    *mapping,
                                    cell_update_flags,
                                    introspection.n_compositional_fields,
                                    stokes_dofs_per_cell,
                                    parameters.include_melt_transport,
                                    rebuild_stokes_matrix),
         internal::Assembly::CopyData::
         StokesPreconditioner<dim> (stokes_dofs_per_cell));

    system_preconditioner_matrix.compress(VectorOperation::add);
  }



  template <int dim>
  void
  Simulator<dim>::build_stokes_preconditioner ()
  {
    if (rebuild_stokes_preconditioner == false)
      return;

    if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_gmg)
      return;
    else if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::block_amg)
      {
        // continue below
      }
    else if (parameters.stokes_solver_type == Parameters<dim>::StokesSolverType::direct_solver)
      return;
    else
      AssertThrow(false, ExcNotImplemented());

    TimerOutput::Scope timer (computing_timer, "Build Stokes preconditioner");
    pcout << "   Rebuilding Stokes preconditioner..." << std::flush;

    // first assemble the raw matrices necessary for the preconditioner
    assemble_stokes_preconditioner ();

    // then extract the other information necessary to build the
    // AMG preconditioners for the A and M blocks
    std::vector<std::vector<bool> > constant_modes;
    DoFTools::extract_constant_modes (dof_handler,
                                      introspection.component_masks.velocities,
                                      constant_modes);

    if (parameters.include_melt_transport)
      Mp_preconditioner = std_cxx14::make_unique<LinearAlgebra::PreconditionAMG>();
    else
      Mp_preconditioner = std_cxx14::make_unique<LinearAlgebra::PreconditionILU>();

    Amg_preconditioner = std_cxx14::make_unique<LinearAlgebra::PreconditionAMG>();

    LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
#ifdef ASPECT_USE_PETSC
    Amg_data.symmetric_operator = false;
#else
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = true;

    // set the AMG parameters in a way that minimizes the run
    // time. compared to some of the deal.II tutorial programs, we
    // found that it pays off to set the aggregation threshold to
    // zero, especially for ill-conditioned problems with large
    // variations in the viscosity
    //
    // for extensive benchmarking of various settings of these
    // parameters and others, see
    // https://github.com/geodynamics/aspect/pull/234
    Amg_data.smoother_type = parameters.AMG_smoother_type.c_str();
    Amg_data.smoother_sweeps = parameters.AMG_smoother_sweeps;
    Amg_data.aggregation_threshold = parameters.AMG_aggregation_threshold;
    Amg_data.output_details = parameters.AMG_output_details;
#endif

    /*  The stabilization term for the free surface (Kaus et. al., 2010)
     *  makes changes to the system matrix which are of the same form as
     *  boundary stresses. If these stresses are not also added to the
     *  system_preconditioner_matrix, then it fails to be very good as a
     *  preconditioner. Instead, we just pass the system_matrix to the
     *  AMG precondition initialization so that it builds the preconditioner
     *  directly from that. However, we still need the mass matrix for the
     *  pressure block which is assembled in the preconditioner matrix.
     *  So rather than build a different preconditioner matrix which only
     *  does the mass matrix, we just reuse the same system_preconditioner_matrix
     *  for the Mp_preconditioner block.  Maybe a bit messy*/

    if (parameters.include_melt_transport == false)
      {
        LinearAlgebra::PreconditionILU *Mp_preconditioner_ILU
          = dynamic_cast<LinearAlgebra::PreconditionILU *> (Mp_preconditioner.get());
        Mp_preconditioner_ILU->initialize (system_preconditioner_matrix.block(1,1));
      }
    else
      {
        // in the case of melt transport we have an AMG preconditioner for the lower right block.
        LinearAlgebra::PreconditionAMG::AdditionalData Amg_data;
#ifdef ASPECT_USE_PETSC
        Amg_data.symmetric_operator = false;
#else
        std::vector<std::vector<bool> > constant_modes;
        dealii::ComponentMask cm_pressure = introspection.component_masks.pressure;
        if (parameters.include_melt_transport)
          cm_pressure = cm_pressure | introspection.variable("compaction pressure").component_mask;
        DoFTools::extract_constant_modes (dof_handler,
                                          cm_pressure,
                                          constant_modes);

        Amg_data.elliptic = true;
        Amg_data.higher_order_elements = false;

        Amg_data.smoother_sweeps = 2;
        Amg_data.coarse_type = "symmetric Gauss-Seidel";
#endif

        LinearAlgebra::PreconditionAMG *Mp_preconditioner_AMG
          = dynamic_cast<LinearAlgebra::PreconditionAMG *> (Mp_preconditioner.get());
        Mp_preconditioner_AMG->initialize (system_preconditioner_matrix.block(1,1), Amg_data);
      }

    if (parameters.use_full_A_block_preconditioner)
      Amg_preconditioner->initialize (system_matrix.block(0,0),
                                      Amg_data);
    else
      Amg_preconditioner->initialize (system_preconditioner_matrix.block(0,0),
                                      Amg_data);

    rebuild_stokes_preconditioner = false;

    pcout << std::endl;
  }



  template <int dim>
  void
  Simulator<dim>::
  local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                internal::Assembly::Scratch::StokesSystem<dim> &scratch,
                                internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    // First get all dof indices of the current cell, then extract those
    // that correspond to the Stokes system we are interested in.
    // Note that assemblers below can modify this list of dofs, if they in fact
    // assemble a different system than the standard Stokes system (e.g. in
    // models with melt transport).

    cell->get_dof_indices (scratch.local_dof_indices);

    data.extract_stokes_dof_indices (scratch.local_dof_indices, introspection, finite_element);

    // Prepare the data structures for assembly
    scratch.reinit(cell);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;
    if (do_pressure_rhs_compatibility_modification)
      data.local_pressure_shape_function_integrals = 0;

    // initialize the material model data on the cell
    const bool update_strain_rate =
      assemble_newton_stokes_system || this->parameters.enable_prescribed_dilation || rebuild_stokes_matrix;
    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         cell,
                                         update_strain_rate,
                                         scratch.material_model_inputs);

    for (unsigned int i=0; i<assemblers->stokes_system.size(); ++i)
      assemblers->stokes_system[i]->create_additional_material_model_outputs(scratch.material_model_outputs);

    material_model->evaluate(scratch.material_model_inputs,
                             scratch.material_model_outputs);
    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.finite_element_values.get_quadrature(),
                                               scratch.finite_element_values.get_mapping(),
                                               scratch.material_model_outputs);

    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.velocity_values);
    if (assemble_newton_stokes_system)
      scratch.finite_element_values[introspection.extractors.velocities].get_function_divergences(current_linearization_point,scratch.velocity_divergence);
    if (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::hydrostatic_compression)
      scratch.finite_element_values[introspection.extractors.temperature].get_function_gradients(current_linearization_point,
          scratch.temperature_gradients);


    const bool use_reference_density_profile = (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::reference_density_profile)
                                               || (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile);
    if (use_reference_density_profile)
      {
        for (unsigned int q=0; q<scratch.finite_element_values.n_quadrature_points; ++q)
          {
            scratch.reference_densities[q] = adiabatic_conditions->density(scratch.material_model_inputs.position[q]);
            scratch.reference_densities_depth_derivative[q] = adiabatic_conditions->density_derivative(scratch.material_model_inputs.position[q]);
          }
      }


    // trigger the invocation of the various functions that actually do
    // all of the assembling
    for (unsigned int i=0; i<assemblers->stokes_system.size(); ++i)
      assemblers->stokes_system[i]->execute(scratch,data);

    if (!assemblers->stokes_system_on_boundary_face.empty())
      {
        // then also work on possible face terms. if necessary, initialize
        // the material model data on faces
        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
          if (cell->at_boundary(face_number) && !cell->has_periodic_neighbor(face_number))
            {
              scratch.reinit(cell, face_number);
              if (assemblers->stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data)
                {
                  const bool need_viscosity = rebuild_stokes_matrix |
                                              assemblers->stokes_system_assembler_on_boundary_face_properties.need_viscosity;

                  compute_material_model_input_values (current_linearization_point,
                                                       scratch.face_finite_element_values,
                                                       cell,
                                                       need_viscosity,
                                                       scratch.face_material_model_inputs);

                  for (unsigned int i=0; i<assemblers->stokes_system_on_boundary_face.size(); ++i)
                    assemblers->stokes_system_on_boundary_face[i]->create_additional_material_model_outputs(scratch.face_material_model_outputs);

                  material_model->evaluate(scratch.face_material_model_inputs,
                                           scratch.face_material_model_outputs);

                  // TODO: Currently we do not supply reference density values to Stokes face assemblers.
                  // This seems acceptable for now, since the only face assemblers are the ones for the melt
                  // assembly where one would not want to use the reference density anyway. In case the reference
                  // density is ever needed, those assemblers can also query the adiabatic conditions for the
                  // reference density.

                  // TODO: the following doesn't currently compile because the get_quadrature() call returns
                  //  a dim-1 dimensional quadrature
                  // MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                  //                                            cell,
                  //                                            compressibility * density *
                  //                                            scratch.face_finite_element_values.get_mapping(),
                  //                                            scratch.face_material_model_outputs);
                }

              for (unsigned int i=0; i<assemblers->stokes_system_on_boundary_face.size(); ++i)
                assemblers->stokes_system_on_boundary_face[i]->execute(scratch,data);
            }
      }
  }



  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_stokes_system (const internal::Assembly::CopyData::StokesSystem<dim> &data)
  {
    if (rebuild_stokes_matrix == true)
      current_constraints.distribute_local_to_global (data.local_matrix,
                                                      data.local_rhs,
                                                      data.local_dof_indices,
                                                      system_matrix,
                                                      system_rhs);
    else
      current_constraints.distribute_local_to_global (data.local_rhs,
                                                      data.local_dof_indices,
                                                      system_rhs);

    if (do_pressure_rhs_compatibility_modification)
      current_constraints.distribute_local_to_global (data.local_pressure_shape_function_integrals,
                                                      data.local_dof_indices,
                                                      pressure_shape_function_integrals);
  }



  template <int dim>
  void Simulator<dim>::assemble_stokes_system ()
  {
    std::string timer_section_name = "Assemble Stokes system";

    if (assemble_newton_stokes_system)
      {
        if (!assemble_newton_stokes_matrix && !stokes_matrix_free)
          timer_section_name += " rhs";
        else if (assemble_newton_stokes_matrix && newton_handler->parameters.newton_derivative_scaling_factor == 0)
          timer_section_name += " Picard";
        else if (assemble_newton_stokes_matrix && newton_handler->parameters.newton_derivative_scaling_factor != 0)
          timer_section_name += " Newton";
      }

    if (stokes_matrix_free)
      {
        rebuild_stokes_matrix = false;
        assemble_newton_stokes_matrix = false;
        timer_section_name += " rhs";
      }

    TimerOutput::Scope timer (computing_timer,
                              timer_section_name);


    if (rebuild_stokes_matrix == true)
      system_matrix = 0;

    // We are using constraints.distribute_local_to_global() without a matrix
    // if we do not rebuild the Stokes matrix. This produces incorrect results
    // when having inhomogeneous constraints. Make sure that we can not have
    // this situation (no active boundary conditions means that only
    // no-slip/free slip are used). This should not happen as we set this up
    // correctly before calling this function.
    // Note that for Dirichlet boundary conditions in matrix-free computations,
    // we will update the right-hand side with boundary information in
    // StokesMatrixFreeHandler::correct_stokes_rhs().
    if (!stokes_matrix_free)
      Assert(rebuild_stokes_matrix || boundary_velocity_manager.get_active_boundary_velocity_conditions().size()==0,
             ExcInternalError("If we have inhomogeneous constraints, we must re-assemble the system matrix."));

    system_rhs = 0;
    if (do_pressure_rhs_compatibility_modification)
      pressure_shape_function_integrals = 0;

    const QGauss<dim>   quadrature_formula(parameters.stokes_velocity_degree+1);
    const QGauss<dim-1> face_quadrature_formula(parameters.stokes_velocity_degree+1);

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    // determine which updates flags we need on cells and faces
    const UpdateFlags cell_update_flags
      = (update_values    |
         update_gradients |
         update_quadrature_points  |
         update_JxW_values)
        |
        assemblers->stokes_system_assembler_properties.needed_update_flags;
    const UpdateFlags face_update_flags
      = (
          // see if we need to assemble traction boundary conditions.
          // only if so do we actually need to have an FEFaceValues object
          parameters.prescribed_traction_boundary_indicators.size() > 0
          ?
          update_values |
          update_quadrature_points |
          update_normal_vectors |
          update_JxW_values
          :
          update_default)
        |
        (assemblers->stokes_system_assembler_on_boundary_face_properties.need_face_material_model_data
         ?
         // if we need a material model input on the faces, we need to
         // also be able to compute the strain rate
         update_gradients
         :
         update_default)
        |
        assemblers->stokes_system_assembler_on_boundary_face_properties.needed_update_flags;

    unsigned int stokes_dofs_per_cell = dim * finite_element.base_element(introspection.base_elements.velocities).dofs_per_cell
                                        + finite_element.base_element(introspection.base_elements.pressure).dofs_per_cell;

    if (parameters.include_melt_transport)
      stokes_dofs_per_cell += finite_element.base_element(introspection.variable("compaction pressure").base_index).dofs_per_cell;

    const bool use_reference_density_profile = (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::reference_density_profile)
                                               || (parameters.formulation_mass_conservation == Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile);

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::Assembly::Scratch::StokesSystem<dim> &scratch,
                      internal::Assembly::CopyData::StokesSystem<dim> &data)
    {
      this->local_assemble_stokes_system(cell, scratch, data);
    };

    auto copier = [&](const internal::Assembly::CopyData::StokesSystem<dim> &data)
    {
      this->copy_local_to_global_stokes_system(data);
    };

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         worker,
         copier,
         internal::Assembly::Scratch::
         StokesSystem<dim> (finite_element, *mapping, quadrature_formula,
                            face_quadrature_formula,
                            cell_update_flags,
                            face_update_flags,
                            introspection.n_compositional_fields,
                            stokes_dofs_per_cell,
                            parameters.include_melt_transport,
                            use_reference_density_profile,
                            rebuild_stokes_matrix,
                            assemble_newton_stokes_matrix),
         internal::Assembly::CopyData::
         StokesSystem<dim> (stokes_dofs_per_cell,
                            do_pressure_rhs_compatibility_modification));

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    // If we change the system_rhs, matrix-free Stokes must update
    if (stokes_matrix_free)
      {
        stokes_matrix_free->evaluate_material_model();
        stokes_matrix_free->correct_stokes_rhs();
      }

    // if the model is compressible then we need to adjust the right hand
    // side of the equation to make it compatible with the matrix on the
    // left
    if (do_pressure_rhs_compatibility_modification)
      {
        pressure_shape_function_integrals.compress(VectorOperation::add);
        make_pressure_rhs_compatible(system_rhs);
      }

    // record that we have just rebuilt the matrix
    rebuild_stokes_matrix = false;
  }



  template <int dim>
  void
  Simulator<dim>::build_advection_preconditioner(const AdvectionField &advection_field,
                                                 LinearAlgebra::PreconditionILU &preconditioner,
                                                 const double diagonal_strengthening)
  {
    TimerOutput::Scope timer (computing_timer, (advection_field.is_temperature() ?
                                                "Build temperature preconditioner" :
                                                "Build composition preconditioner"));

    const unsigned int block_idx = advection_field.block_index(introspection);


    LinearAlgebra::PreconditionILU::AdditionalData data;
    data.ilu_atol = diagonal_strengthening;

    preconditioner.initialize (system_matrix.block(block_idx, block_idx), data);
  }



  template <int dim>
  void Simulator<dim>::
  local_assemble_advection_system (const AdvectionField     &advection_field,
                                   const Vector<double>           &viscosity_per_cell,
                                   const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                                   internal::Assembly::CopyData::AdvectionSystem<dim> &data)
  {
    // also have the number of dofs that correspond just to the element for
    // the system we are currently trying to assemble
    const unsigned int advection_dofs_per_cell = data.local_dof_indices.size();

    Assert (advection_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
    Assert (scratch.grad_phi_field.size() == advection_dofs_per_cell, ExcInternalError());
    Assert (scratch.phi_field.size() == advection_dofs_per_cell, ExcInternalError());

    const FEValuesExtractors::Scalar solution_field = advection_field.scalar_extractor(introspection);

    const unsigned int solution_component = advection_field.component_index(introspection);

    scratch.reinit(cell);

    // get all dof indices on the current cell, then extract those
    // that correspond to the solution_field we are interested in
    cell->get_dof_indices (scratch.local_dof_indices);
    for (unsigned int i=0, i_advection=0; i_advection<advection_dofs_per_cell; /*increment at end of loop*/)
      {
        if (finite_element.system_to_component_index(i).first == solution_component)
          {
            data.local_dof_indices[i_advection] = scratch.local_dof_indices[i];
            ++i_advection;
          }
        ++i;
      }

    data.local_matrix = 0;
    data.local_rhs = 0;

    scratch.finite_element_values[solution_field].get_function_values (old_solution,
                                                                       scratch.old_field_values);
    scratch.finite_element_values[solution_field].get_function_values (old_old_solution,
                                                                       scratch.old_old_field_values);


    scratch.finite_element_values[introspection.extractors.velocities].get_function_values(current_linearization_point,
        scratch.current_velocity_values);

    if (parameters.include_melt_transport)
      scratch.finite_element_values[introspection.extractors.velocities].get_function_divergences(current_linearization_point,
          scratch.current_velocity_divergences);

    // get the mesh velocity, as we need to subtract it off of the advection systems
    if (parameters.mesh_deformation_enabled)
      scratch.finite_element_values[introspection.extractors.velocities].get_function_values(mesh_deformation->mesh_velocity,
          scratch.mesh_velocity_values);

    // compute material properties and heating terms
    compute_material_model_input_values (current_linearization_point,
                                         scratch.finite_element_values,
                                         cell,
                                         true,
                                         scratch.material_model_inputs);

    for (unsigned int i=0; i<assemblers->advection_system.size(); ++i)
      assemblers->advection_system[i]->create_additional_material_model_outputs(scratch.material_model_outputs);

    heating_model_manager.create_additional_material_model_inputs_and_outputs(scratch.material_model_inputs,
                                                                              scratch.material_model_outputs);
    material_model->fill_additional_material_model_inputs(scratch.material_model_inputs,
                                                          current_linearization_point,
                                                          scratch.finite_element_values,
                                                          introspection);

    material_model->evaluate(scratch.material_model_inputs,
                             scratch.material_model_outputs);
    if (parameters.formulation_temperature_equation ==
        Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
      {
        const unsigned int n_q_points = scratch.finite_element_values.n_quadrature_points;
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            scratch.material_model_outputs.densities[q] = adiabatic_conditions->density(scratch.material_model_inputs.position[q]);
          }
      }

#ifdef DEBUG
    // make sure that if the model does not use operator splitting,
    // the material model outputs do not fill the reaction_rates (because the reaction_terms are used instead)
    if (!parameters.use_operator_splitting)
      {
        material_model->create_additional_named_outputs(scratch.material_model_outputs);
        MaterialModel::ReactionRateOutputs<dim> *reaction_rate_outputs
          = scratch.material_model_outputs.template get_additional_output<MaterialModel::ReactionRateOutputs<dim> >();

        Assert(reaction_rate_outputs == nullptr,
               ExcMessage("You are using a material model where the reaction rate outputs "
                          "are created even though the operator splitting solver option is "
                          "not used in the model, this is not supported! "
                          "If operator splitting is disabled, the reaction_rates should not "
                          "be created at all. If you want to run a model where reactions are "
                          "much faster than the advection, which is what the reaction rate "
                          "outputs are designed for, you should enable operator splitting."));
      }
#endif

    MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                                               cell,
                                               scratch.finite_element_values.get_quadrature(),
                                               scratch.finite_element_values.get_mapping(),
                                               scratch.material_model_outputs);

    heating_model_manager.evaluate(scratch.material_model_inputs,
                                   scratch.material_model_outputs,
                                   scratch.heating_model_outputs);

    // TODO: Compute artificial viscosity once per timestep instead of each time
    // temperature system is assembled (as this might happen more than once per
    // timestep for iterative solvers)
    scratch.artificial_viscosity = viscosity_per_cell[cell->active_cell_index()];
    Assert (scratch.artificial_viscosity >= 0, ExcMessage ("The artificial viscosity needs to be a non-negative quantity."));

    // trigger the invocation of the various functions that actually do
    // all of the assembling
    for (unsigned int i=0; i<assemblers->advection_system.size(); ++i)
      assemblers->advection_system[i]->execute(scratch,data);

    // then also work on possible face terms. if necessary, initialize
    // the material model data on faces
    const bool has_boundary_face_assemblers = !assemblers->advection_system_on_boundary_face.empty()
                                              && assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;
    const bool has_interior_face_assemblers = !assemblers->advection_system_on_interior_face.empty()
                                              && assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;

    // skip the remainder if no work needs to be done on faces
    if (!has_boundary_face_assemblers && !has_interior_face_assemblers)
      return;

    if (has_interior_face_assemblers)
      {
        // for interior face contributions loop over all possible
        // subfaces of the cell, and reset their matrices.
        for (unsigned int f = 0; f < GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell; ++f)
          {
            data.local_matrices_ext_int[f] = 0;
            data.local_matrices_int_ext[f] = 0;
            data.local_matrices_ext_ext[f] = 0;
            data.assembled_matrices[f] = false;
          }
      }

    for (scratch.face_number=0; scratch.face_number<GeometryInfo<dim>::faces_per_cell; ++scratch.face_number)
      {
        typename DoFHandler<dim>::face_iterator face = cell->face (scratch.face_number);

        if ((has_boundary_face_assemblers && face->at_boundary()) ||
            (has_interior_face_assemblers && !face->at_boundary()))
          {
            (*scratch.face_finite_element_values).reinit (cell, scratch.face_number);

            (*scratch.face_finite_element_values)[introspection.extractors.velocities].get_function_values(current_linearization_point,
                scratch.face_current_velocity_values);

            // get the mesh velocity, as we need to subtract it off of the advection systems
            if (parameters.mesh_deformation_enabled)
              (*scratch.face_finite_element_values)[introspection.extractors.velocities].get_function_values(mesh_deformation->mesh_velocity,
                  scratch.face_mesh_velocity_values);

            if (assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_material_model_data)
              {
                compute_material_model_input_values (current_linearization_point,
                                                     *scratch.face_finite_element_values,
                                                     cell,
                                                     true,
                                                     scratch.face_material_model_inputs);

                for (unsigned int i=0; i<assemblers->advection_system_on_boundary_face.size(); ++i)
                  assemblers->advection_system_on_boundary_face[i]->create_additional_material_model_outputs(scratch.face_material_model_outputs);

                for (unsigned int i=0; i<assemblers->advection_system_on_interior_face.size(); ++i)
                  assemblers->advection_system_on_interior_face[i]->create_additional_material_model_outputs(scratch.face_material_model_outputs);

                heating_model_manager.create_additional_material_model_inputs_and_outputs(scratch.face_material_model_inputs,
                                                                                          scratch.face_material_model_outputs);
                material_model->fill_additional_material_model_inputs(scratch.face_material_model_inputs,
                                                                      current_linearization_point,
                                                                      *scratch.face_finite_element_values,
                                                                      introspection);

                material_model->evaluate(scratch.face_material_model_inputs,
                                         scratch.face_material_model_outputs);

                if (parameters.formulation_temperature_equation ==
                    Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                  {
                    const unsigned int n_q_points = scratch.face_finite_element_values->n_quadrature_points;
                    for (unsigned int q=0; q<n_q_points; ++q)
                      {
                        scratch.face_material_model_outputs.densities[q] = adiabatic_conditions->density(scratch.face_material_model_inputs.position[q]);
                      }
                  }

                heating_model_manager.evaluate(scratch.face_material_model_inputs,
                                               scratch.face_material_model_outputs,
                                               scratch.face_heating_model_outputs);
                // TODO: the following doesn't currently compile because the get_quadrature() call returns
                //  a dim-1 dimensional quadrature
                // MaterialModel::MaterialAveraging::average (parameters.material_averaging,
                //                                            cell,
                //                                            scratch.face_finite_element_values.get_quadrature(),
                //                                            scratch.face_finite_element_values.get_mapping(),
                //                                            scratch.face_material_model_outputs);
              }

            // handle faces at periodic boundaries like interior faces
            if (face->at_boundary() && !cell->has_periodic_neighbor (scratch.face_number))
              for (unsigned int i=0; i<assemblers->advection_system_on_boundary_face.size(); ++i)
                assemblers->advection_system_on_boundary_face[i]->execute(scratch,data);
            else
              for (unsigned int i=0; i<assemblers->advection_system_on_interior_face.size(); ++i)
                assemblers->advection_system_on_interior_face[i]->execute(scratch,data);
          }
      }
  }

  template <int dim>
  void
  Simulator<dim>::
  copy_local_to_global_advection_system (const AdvectionField &advection_field,
                                         const internal::Assembly::CopyData::AdvectionSystem<dim> &data)
  {
    // copy entries into the global matrix. note that these local contributions
    // only correspond to the advection dofs, as assembled above
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);

    /* In the following, we copy DG contributions element by element. This
     * is allowed since there are no constraints imposed on discontinuous fields.
     */
    if (!assemblers->advection_system_on_interior_face.empty() &&
        assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation)
      {
        for (unsigned int f=0; f<GeometryInfo<dim>::max_children_per_face
             * GeometryInfo<dim>::faces_per_cell; ++f)
          {
            if (data.assembled_matrices[f])
              {
                for (unsigned int i=0; i<data.local_dof_indices.size(); ++i)
                  for (unsigned int j=0; j<data.neighbor_dof_indices[f].size(); ++j)
                    {
                      system_matrix.add (data.local_dof_indices[i],
                                         data.neighbor_dof_indices[f][j],
                                         data.local_matrices_int_ext[f](i,j));
                      system_matrix.add (data.neighbor_dof_indices[f][j],
                                         data.local_dof_indices[i],
                                         data.local_matrices_ext_int[f](j,i));
                    }

                for (unsigned int i=0; i<data.neighbor_dof_indices[f].size(); ++i)
                  for (unsigned int j=0; j<data.neighbor_dof_indices[f].size(); ++j)
                    system_matrix.add (data.neighbor_dof_indices[f][i],
                                       data.neighbor_dof_indices[f][j],
                                       data.local_matrices_ext_ext[f](i,j));
              }
          }
      }
  }



  template <int dim>
  void Simulator<dim>::assemble_advection_system (const AdvectionField &advection_field)
  {
    TimerOutput::Scope timer (computing_timer, (advection_field.is_temperature() ?
                                                "Assemble temperature system" :
                                                "Assemble composition system"));

    const unsigned int block_idx = advection_field.block_index(introspection);

    if (!advection_field.is_temperature() && advection_field.compositional_variable!=0)
      {
        // Allocate the system matrix for the current compositional field by
        // reusing the Trilinos sparsity pattern from the matrix stored for
        // composition 0 (this is the place we allocate the matrix at).
        const unsigned int block0_idx = AdvectionField::composition(0).block_index(introspection);
        system_matrix.block(block_idx, block_idx).reinit(system_matrix.block(block0_idx, block0_idx));
      }

    system_matrix.block(block_idx, block_idx) = 0;
    system_rhs.block(block_idx) = 0;


    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

    Vector<double> viscosity_per_cell;
    viscosity_per_cell.reinit(triangulation.n_active_cells());
    get_artificial_viscosity(viscosity_per_cell, advection_field);

    // We have to assemble the term u.grad phi_i * phi_j, which is
    // of total polynomial degree
    //   stokes_deg + 2*temp_deg -1
    // (or similar for comp_deg). This suggests using a Gauss
    // quadrature formula of order
    //   temp_deg + stokes_deg/2
    // rounded up (note that x/2 rounded up
    // equals (x+1)/2 using integer division.)
    //
    // (Note: All compositional fields have the same base element and therefore
    // the same composition_degree. Thus, we do not need to find out the degree
    // of the current field, but use the global instead)
    const unsigned int advection_quadrature_degree = advection_field.polynomial_degree(introspection)
                                                     +
                                                     (parameters.stokes_velocity_degree+1)/2;

    const bool allocate_face_quadrature = (!assemblers->advection_system_on_boundary_face.empty() ||
                                           !assemblers->advection_system_on_interior_face.empty()) &&
                                          assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;
    const bool allocate_neighbor_contributions = !assemblers->advection_system_on_interior_face.empty() &&
                                                 assemblers->advection_system_assembler_on_face_properties[advection_field.field_index()].need_face_finite_element_evaluation;;

    const bool use_supg = (parameters.advection_stabilization_method
                           == Parameters<dim>::AdvectionStabilizationMethod::supg);

    // When using SUPG, we need to compute hessians to be able to compute the residual:
    const UpdateFlags update_flags = update_values |
                                     update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values |
                                     ((use_supg) ? update_hessians : UpdateFlags(0));

    const UpdateFlags face_update_flags = (allocate_face_quadrature ?
                                           update_values |
                                           update_gradients |
                                           update_quadrature_points |
                                           update_normal_vectors |
                                           update_JxW_values
                                           :
                                           update_default);

    auto worker = [&](const typename DoFHandler<dim>::active_cell_iterator &cell,
                      internal::Assembly::Scratch::AdvectionSystem<dim> &scratch,
                      internal::Assembly::CopyData::AdvectionSystem<dim> &data)
    {
      this->local_assemble_advection_system(advection_field, viscosity_per_cell, cell, scratch, data);
    };

    auto copier = [&](const internal::Assembly::CopyData::AdvectionSystem<dim> &data)
    {
      this->copy_local_to_global_advection_system(advection_field, data);
    };

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         worker,
         copier,
         internal::Assembly::Scratch::
         AdvectionSystem<dim> (finite_element,
                               finite_element.base_element(advection_field.base_element(introspection)),
                               *mapping,
                               QGauss<dim>(advection_quadrature_degree),
                               /* Only generate a valid face quadrature if necessary.
                                * Otherwise, generate invalid face quadrature rule.
                                */
                               (allocate_face_quadrature ?
                                QGauss<dim-1>(advection_quadrature_degree) :
                                Quadrature<dim-1> ()),
                               update_flags,
                               face_update_flags,
                               introspection.n_compositional_fields,
                               advection_field),
         internal::Assembly::CopyData::
         AdvectionSystem<dim> (finite_element.base_element(advection_field.base_element(introspection)),
                               allocate_neighbor_contributions));

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }
}



// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::set_assemblers (); \
  template void Simulator<dim>::local_assemble_stokes_preconditioner ( \
                                                                       const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                       internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch, \
                                                                       internal::Assembly::CopyData::StokesPreconditioner<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_stokes_preconditioner ( \
                                                                             const internal::Assembly::CopyData::StokesPreconditioner<dim> &data); \
  template void Simulator<dim>::assemble_stokes_preconditioner (); \
  template void Simulator<dim>::build_stokes_preconditioner (); \
  template void Simulator<dim>::local_assemble_stokes_system ( \
                                                               const DoFHandler<dim>::active_cell_iterator &cell, \
                                                               internal::Assembly::Scratch::StokesSystem<dim>  &scratch, \
                                                               internal::Assembly::CopyData::StokesSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_stokes_system ( \
                                                                     const internal::Assembly::CopyData::StokesSystem<dim> &data); \
  template void Simulator<dim>::assemble_stokes_system (); \
  template void Simulator<dim>::build_advection_preconditioner (const AdvectionField &, \
                                                                aspect::LinearAlgebra::PreconditionILU &preconditioner, \
                                                                const double diagonal_strengthening); \
  template void Simulator<dim>::local_assemble_advection_system ( \
                                                                  const AdvectionField          &advection_field, \
                                                                  const Vector<double>           &viscosity_per_cell, \
                                                                  const DoFHandler<dim>::active_cell_iterator &cell, \
                                                                  internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch, \
                                                                  internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_advection_system ( \
                                                                        const AdvectionField          &advection_field, \
                                                                        const internal::Assembly::CopyData::AdvectionSystem<dim> &data); \
  template void Simulator<dim>::assemble_advection_system (const AdvectionField     &advection_field); \
  template void Simulator<dim>::compute_material_model_input_values ( \
                                                                      const LinearAlgebra::BlockVector                      &input_solution, \
                                                                      const FEValuesBase<dim,dim>                           &input_finite_element_values, \
                                                                      const DoFHandler<dim>::active_cell_iterator  &cell, \
                                                                      const bool                                             compute_strainrate, \
                                                                      MaterialModel::MaterialModelInputs<dim>               &material_model_inputs) const;



  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
