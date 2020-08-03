/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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
#include <aspect/parameters.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/mesh_refinement/volume_of_fluid_interface.h>

#include <deal.II/base/work_stream.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_cartesian.h>

#ifdef ASPECT_USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

using namespace dealii;

namespace aspect
{
  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        VolumeOfFluidSystem<dim>::VolumeOfFluidSystem (const FiniteElement<dim> &finite_element,
                                                       const FiniteElement<dim> &volume_of_fluid_element,
                                                       const Mapping<dim>       &mapping,
                                                       const Quadrature<dim>    &quadrature,
                                                       const Quadrature<dim-1>  &face_quadrature)
          :
          finite_element_values (mapping,
                                 finite_element, quadrature,
                                 update_values |
                                 update_gradients |
                                 update_JxW_values),
          neighbor_finite_element_values (mapping,
                                          finite_element, quadrature,
                                          update_values |
                                          update_gradients |
                                          update_JxW_values),
          face_finite_element_values (mapping,
                                      finite_element, face_quadrature,
                                      update_values |
                                      update_quadrature_points |
                                      update_gradients |
                                      update_normal_vectors |
                                      update_JxW_values),
          neighbor_face_finite_element_values (mapping,
                                               finite_element, face_quadrature,
                                               update_values |
                                               update_gradients |
                                               update_normal_vectors |
                                               update_JxW_values),
          subface_finite_element_values (mapping,
                                         finite_element, face_quadrature,
                                         update_values |
                                         update_gradients |
                                         update_normal_vectors |
                                         update_JxW_values),
          local_dof_indices(finite_element.dofs_per_cell),
          phi_field (volume_of_fluid_element.dofs_per_cell, numbers::signaling_nan<double>()),
          old_field_values (quadrature.size(), numbers::signaling_nan<double>()),
          cell_i_n_values (quadrature.size(), numbers::signaling_nan<Tensor<1, dim> > ()),
          cell_i_d_values (quadrature.size(), numbers::signaling_nan<double> ()),
          face_current_velocity_values (face_quadrature.size(), numbers::signaling_nan<Tensor<1, dim> >()),
          face_old_velocity_values (face_quadrature.size(), numbers::signaling_nan<Tensor<1, dim> >()),
          face_old_old_velocity_values (face_quadrature.size(), numbers::signaling_nan<Tensor<1, dim> >()),
          neighbor_old_values (face_quadrature.size(), numbers::signaling_nan<double>()),
          neighbor_i_n_values (face_quadrature.size(), numbers::signaling_nan<Tensor<1, dim> >()),
          neighbor_i_d_values (face_quadrature.size(), numbers::signaling_nan<double>())
        {}



        template <int dim>
        VolumeOfFluidSystem<dim>::VolumeOfFluidSystem (const VolumeOfFluidSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),
          neighbor_finite_element_values (scratch.neighbor_finite_element_values.get_mapping(),
                                          scratch.neighbor_finite_element_values.get_fe(),
                                          scratch.neighbor_finite_element_values.get_quadrature(),
                                          scratch.neighbor_finite_element_values.get_update_flags()),
          face_finite_element_values (scratch.face_finite_element_values.get_mapping(),
                                      scratch.face_finite_element_values.get_fe(),
                                      scratch.face_finite_element_values.get_quadrature(),
                                      scratch.face_finite_element_values.get_update_flags()),
          neighbor_face_finite_element_values (scratch.neighbor_face_finite_element_values.get_mapping(),
                                               scratch.neighbor_face_finite_element_values.get_fe(),
                                               scratch.neighbor_face_finite_element_values.get_quadrature(),
                                               scratch.neighbor_face_finite_element_values.get_update_flags()),
          subface_finite_element_values (scratch.subface_finite_element_values.get_mapping(),
                                         scratch.subface_finite_element_values.get_fe(),
                                         scratch.subface_finite_element_values.get_quadrature(),
                                         scratch.subface_finite_element_values.get_update_flags()),
          local_dof_indices (scratch.finite_element_values.get_fe().dofs_per_cell),
          phi_field (scratch.phi_field),
          old_field_values (scratch.old_field_values),
          cell_i_n_values (scratch.cell_i_n_values),
          cell_i_d_values (scratch.cell_i_d_values),
          face_current_velocity_values (scratch.face_current_velocity_values),
          face_old_velocity_values (scratch.face_old_velocity_values),
          face_old_old_velocity_values (scratch.face_old_old_velocity_values),
          neighbor_old_values (scratch.neighbor_old_values),
          neighbor_i_n_values (scratch.neighbor_i_n_values),
          neighbor_i_d_values (scratch.neighbor_i_d_values)
        {}
      }



      namespace CopyData
      {
        template <int dim>
        VolumeOfFluidSystem<dim>::VolumeOfFluidSystem(const FiniteElement<dim> &finite_element)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_rhs (finite_element.dofs_per_cell),
          local_dof_indices (finite_element.dofs_per_cell)
        {
          TableIndices<2> mat_size(finite_element.dofs_per_cell,
                                   finite_element.dofs_per_cell);
          for (unsigned int i=0;
               i < GeometryInfo<dim>::max_children_per_face *GeometryInfo<dim>::faces_per_cell;
               ++i)
            {
              face_contributions_mask[i] = false;
              local_face_rhs[i].reinit (finite_element.dofs_per_cell);
              local_face_matrices_ext_ext[i].reinit(mat_size);
              neighbor_dof_indices[i].resize(finite_element.dofs_per_cell);
            }
        }



        template<int dim>
        VolumeOfFluidSystem<dim>::VolumeOfFluidSystem(const VolumeOfFluidSystem &data)
          :
          local_matrix (data.local_matrix),
          local_rhs (data.local_rhs),
          local_face_rhs (data.local_face_rhs),
          local_face_matrices_ext_ext (data.local_face_matrices_ext_ext),
          local_dof_indices (data.local_dof_indices),
          neighbor_dof_indices (data.neighbor_dof_indices)
        {
          for (unsigned int i=0;
               i < GeometryInfo<dim>::max_children_per_face *GeometryInfo<dim>::faces_per_cell;
               ++i)
            {
              face_contributions_mask[i] = false;
            }
        }
      }
    }
  }



  template <int dim>
  VolumeOfFluidField<dim>::VolumeOfFluidField(const FEVariable<dim> &volume_fraction,
                                              const FEVariable<dim> &reconstruction,
                                              const FEVariable<dim> &level_set,
                                              const unsigned int composition_index)
    : volume_fraction (volume_fraction),
      reconstruction (reconstruction),
      level_set (level_set),
      composition_index(composition_index)
  {}



  template <int dim>
  VolumeOfFluidHandler<dim>::VolumeOfFluidHandler (Simulator<dim> &simulator,
                                                   ParameterHandler &prm)
    : sim (simulator),
      assembler ()
  {
    this->initialize_simulator(sim);
    assembler.initialize_simulator(sim);
    parse_parameters (prm);

    this->get_signals().edit_finite_element_variables.connect(std::bind(&aspect::VolumeOfFluidHandler<dim>::edit_finite_element_variables,
                                                                        std::ref(*this),
                                                                        std::placeholders::_1));
    this->get_signals().post_set_initial_state.connect(std::bind(&aspect::VolumeOfFluidHandler<dim>::set_initial_volume_fractions,
                                                                 std::ref(*this)));
  }



  template <int dim>
  void
  VolumeOfFluidHandler<dim>::edit_finite_element_variables (std::vector<VariableDeclaration<dim> > &vars)
  {
    for (unsigned int f=0; f<n_volume_of_fluid_fields; ++f)
      {
        // Add declaration for volume fraction field
        vars.push_back(VariableDeclaration<dim>("volume_fraction_"+volume_of_fluid_field_names[f],
                                                std::unique_ptr<FiniteElement<dim>>(
                                                  new FE_DGQ<dim>(0)),
                                                1,
                                                1));

        // Add declaration for reconstructed interface cache
        vars.push_back(VariableDeclaration<dim>("volume_of_fluid_interface_reconstruction_"+volume_of_fluid_field_names[f],
                                                std::unique_ptr<FiniteElement<dim>>(
                                                  new FE_DGQ<dim>(0)),
                                                dim+1,
                                                1));

        vars.push_back(VariableDeclaration<dim>("volume_of_fluid_contour_"+volume_of_fluid_field_names[f],
                                                std::unique_ptr<FiniteElement<dim>>(
                                                  new FE_DGQ<dim>(1)),
                                                1,
                                                1));
      }
  }



  template <int dim>
  void
  VolumeOfFluidHandler<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Volume of Fluid");
    {
      prm.declare_entry ("Volume fraction threshold", "1e-6",
                         Patterns::Double (0., 1.),
                         "Minimum significant volume. Fluid fractions below this value are considered to be zero.");

      prm.declare_entry ("Volume of Fluid solver tolerance", "1e-12",
                         Patterns::Double(0., 1.),
                         "The relative tolerance up to which the linear system "
                         "for the Volume of Fluid system gets solved. See "
                         "'Solver parameters/Composition solver tolerance' "
                         "for more details.");

      prm.declare_entry ("Number initialization samples", "3",
                         Patterns::Integer (1),
                         "Number of divisions per dimension when computing the initial volume fractions."
                         "If set to the default of 3 for a 2D model, then initialization will be based on "
                         "the initialization criterion at $3^2=9$ points within each cell. If the initialization "
                         "based on a composition style initial condition, a larger value may be desired for better "
                         "approximation of the initial fluid fractions. Smaller values will suffice in the case of "
                         "level set initializations due to the presence of more information to better approximate "
                         "the initial fluid fractions.");
    }
    prm.leave_subsection ();

    prm.enter_subsection("Initial composition model");
    {
      prm.declare_entry("Volume of fluid intialization type", "",
                        Patterns::Map (Patterns::Anything(),
                                       Patterns::Selection("composition|level set")),
                        "A comma separated list denoting the method to be used to "
                        "initialize a composition field specified to be advected using "
                        "the volume of fluid method.\n\n"
                        "The format of valid entries for this parameter is that "
                        "of a map given as ``key1:value1, key2:value2`` where "
                        "each key must be the name of a compositional field "
                        "using the volume of fluid advection method, and the "
                        "value is one of ``composition`` or ``level "
                        "set``. ``composition`` is the default\n\n"
                        "When ``composition is specified, the initial model is "
                        "treated as a standard composition field with bounds "
                        "between 0 and 1 assumed, The initial fluid fractions "
                        "are then based on an iterated midpoint quadrature. "
                        "Resultant volume fractions outside of the bounds will be "
                        "coerced to the nearest valid value (ie 0 or 1). "
                        "If ``level set`` is specified, the intial data will be assumed to "
                        "be in the form of a signed distance level set function "
                        "(i.e. a function which is positive when in the "
                        "fluid, negative outside, and zero on the interface "
                        "and the magnitude is always the distance to the "
                        "interface so the gradient is one everywhere).");
    }
    prm.leave_subsection();
  }



  template <int dim>
  void
  VolumeOfFluidHandler<dim>::parse_parameters (ParameterHandler &prm)
  {
    // Get parameter data
    n_volume_of_fluid_fields = 0;
    std::vector<std::string> names_of_compositional_fields = this->get_parameters().names_of_compositional_fields;
    std::vector<typename Parameters<dim>::AdvectionFieldMethod::Kind> compositional_field_methods = this->get_parameters().compositional_field_methods;

    for (unsigned int i=0; i<names_of_compositional_fields.size(); ++i)
      {
        if (compositional_field_methods[i] == Parameters<dim>::AdvectionFieldMethod::volume_of_fluid)
          {
            // Add this field as the next volume of fluid field
            volume_of_fluid_field_names.push_back(names_of_compositional_fields[i]);
            // Note that compositional field indicies include temperature as field 0, so increase index by 1
            volume_of_fluid_composition_map_index[i+1] = n_volume_of_fluid_fields;
            ++n_volume_of_fluid_fields;
          }
      }

    prm.enter_subsection ("Volume of Fluid");
    {
      volume_fraction_threshold = prm.get_double("Volume fraction threshold");

      volume_of_fluid_solver_tolerance = prm.get_double("Volume of Fluid solver tolerance");

      n_init_samples = prm.get_integer ("Number initialization samples");
    }
    prm.leave_subsection ();

    prm.enter_subsection("Initial composition model");
    {
      const std::vector<std::string> x_initialization_type = Utilities::split_string_list(prm.get("Volume of fluid intialization type"));

      initialization_data_type = std::vector<VolumeOfFluid::VolumeOfFluidInputType::Kind> (n_volume_of_fluid_fields,
                                 VolumeOfFluid::VolumeOfFluidInputType::composition);

      for (std::vector<std::string>::const_iterator p = x_initialization_type.begin();
           p != x_initialization_type.end(); ++p)
        {
          // each entry has the format (white space is optional):
          // <name> : <value (might have spaces)>
          //
          // first tease apart the two halves
          const std::vector<std::string> split_parts = Utilities::split_string_list (*p, ':');
          AssertThrow (split_parts.size() == 2,
                       ExcMessage ("The format for "
                                   "<Initial composition model/Volume of Fluid intialization method> "
                                   "volume of fluid initialization met "
                                   "requires that each entry " "has the form "
                                   "`<name of field> : <method>', but this "
                                   "does not match the number of colons in the "
                                   "entry <" + *p + ">."));

          // get the name of the compositional field
          const std::string key = split_parts[0];

          // check that the names used are actually names of fields,
          // are solved by volume of fluid fields, and are unique in this list
          std::vector<std::string>::iterator field_name_iterator = std::find(names_of_compositional_fields.begin(),
                                                                             names_of_compositional_fields.end(), key);
          AssertThrow (field_name_iterator
                       != names_of_compositional_fields.end(),
                       ExcMessage ("Name of field <" + key +
                                   "> appears in the parameter "
                                   "<Initial composition model/Volume of Fluid intialization method>, but "
                                   "there is no field with this name."));

          const unsigned int compositional_field_index = std::distance(names_of_compositional_fields.begin(),
                                                                       field_name_iterator);

          AssertThrow (compositional_field_methods[compositional_field_index]
                       == Parameters<dim>::AdvectionFieldMethod::volume_of_fluid,
                       ExcMessage ("The field <" + key + "> appears in the parameter "
                                   "<Initial composition model/Volume of Fluid intialization method>, "
                                   "but is not advected by a particle method."));

          AssertThrow (std::count(names_of_compositional_fields.begin(),
                                  names_of_compositional_fields.end(), key) == 1,
                       ExcMessage ("Name of field <" + key + "> appears more "
                                   "than once in the parameter "
                                   "<Initial composition model/Volume of Fluid intialization type>."));

          // Get specification for how to treat initializing data
          const std::string value = split_parts[1];

          if (value == "composition")
            initialization_data_type[volume_of_fluid_composition_map_index[compositional_field_index+1]]
              = VolumeOfFluid::VolumeOfFluidInputType::composition;
          else if (value == "level set")
            initialization_data_type[volume_of_fluid_composition_map_index[compositional_field_index+1]]
              = VolumeOfFluid::VolumeOfFluidInputType::level_set;
          else
            AssertThrow(false,ExcNotImplemented());
        }
    }
    prm.leave_subsection();
  }



  template <int dim>
  void
  VolumeOfFluidHandler<dim>::initialize (ParameterHandler &/*prm*/)
  {
    // Do checks on required assumptions
    AssertThrow(dim==2,
                ExcMessage("Volume of Fluid Interface Tracking is currently only functional for dim=2."));

    AssertThrow(this->get_parameters().CFL_number < 1.0,
                ExcMessage("Volume of Fluid Interface Tracking requires CFL < 1."));

    AssertThrow(!this->get_material_model().is_compressible(),
                ExcMessage("Volume of Fluid Interface Tracking currently assumes incompressiblity."));

    AssertThrow(dynamic_cast<const MappingCartesian<dim> *>(&(this->get_mapping())),
                ExcMessage("Volume of Fluid Interface Tracking currently requires Cartesian Mappings"));

    AssertThrow(!this->get_parameters().mesh_deformation_enabled,
                ExcMessage("Volume of Fluid Interface Tracking is currently incompatible with the Free Surface implementation."));

    AssertThrow(!this->get_parameters().include_melt_transport,
                ExcMessage("Volume of Fluid Interface Tracking has not been tested with melt transport yet, so inclusion of both is currently disabled."))

    if ( this->get_parameters().initial_adaptive_refinement > 0 ||
         this->get_parameters().adaptive_refinement_interval > 0 )
      {
        AssertThrow(this->get_mesh_refinement_manager().template has_matching_mesh_refinement_strategy<MeshRefinement::VolumeOfFluidInterface<dim> >(),
                    ExcMessage("Volume of Fluid Interface Tracking requires that the 'volume of fluid interface' strategy be used for AMR"));

        AssertThrow(this->get_parameters().adaptive_refinement_interval <(1/this->get_parameters().CFL_number),
                    ExcMessage("Currently, the Volume of Fluid Interface Tracking requires that the AMR interval be less than "
                               "1.0/CFL_number to preserve the intended accuracy. Circumventing this limitation requires "
                               "additional care in the reconstruction algorithm."));
      }

    // Gather the created volume fraction data into a structure for easier programmatic referencing

    for (unsigned int f=0; f<n_volume_of_fluid_fields; ++f)
      {
        data.push_back(VolumeOfFluidField<dim>(this->introspection().variable("volume_fraction_"+volume_of_fluid_field_names[f]),
                                               this->introspection().variable("volume_of_fluid_interface_reconstruction_"+volume_of_fluid_field_names[f]),
                                               this->introspection().variable("volume_of_fluid_contour_"+volume_of_fluid_field_names[f]),
                                               this->introspection().compositional_index_for_name(volume_of_fluid_field_names[f])));
      }
  }



  template <int dim>
  unsigned int VolumeOfFluidHandler<dim>::get_n_fields() const
  {
    return n_volume_of_fluid_fields;
  }



  template <int dim>
  const std::string VolumeOfFluidHandler<dim>::name_for_field_index(unsigned int field) const
  {
    Assert(field < n_volume_of_fluid_fields,
           ExcMessage("Invalid field index"));
    return volume_of_fluid_field_names[field];
  }



  template <int dim>
  double VolumeOfFluidHandler<dim>::get_volume_fraction_threshold() const
  {
    return volume_fraction_threshold;
  }



  template <int dim>
  const VolumeOfFluidField<dim> &VolumeOfFluidHandler<dim>::field_struct_for_field_index(unsigned int field) const
  {
    Assert(field < n_volume_of_fluid_fields,
           ExcMessage("Invalid field index"));
    return data[field];
  }



  template <int dim>
  unsigned int VolumeOfFluidHandler<dim>::field_index_for_name(const std::string &composition_fieldname) const
  {
    const unsigned int composition_index = this->introspection().compositional_index_for_name(composition_fieldname);
    if (volume_of_fluid_composition_map_index.count(composition_index) ==0)
      return n_volume_of_fluid_fields;
    return volume_of_fluid_composition_map_index.at(composition_index);
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::do_volume_of_fluid_update (const typename Simulator<dim>::AdvectionField &advection_field)
  {
    const bool direction_order_descending = (this->get_timestep_number() % 2) == 1;
    const VolumeOfFluidField<dim> volume_of_fluid_field = data[volume_of_fluid_composition_map_index[advection_field.field_index()]];

    const unsigned int volume_of_fluid_block_idx = volume_of_fluid_field.volume_fraction.block_index;
    const unsigned int volume_of_fluidN_block_idx = volume_of_fluid_field.reconstruction.block_index;

    // Due to dimensionally split formulation, use Strang (second-order dimensional) splitting
    for (unsigned int direction = 0; direction < dim; ++direction)
      {
        // Only reference old_solution for data from prior substep if this is the first
        // substep for dimensional splitting
        bool update_from_old = (direction == 0);
        // Update base to intermediate solution
        if (!direction_order_descending)
          {
            assemble_volume_of_fluid_system(volume_of_fluid_field, direction, update_from_old);
          }
        else
          {
            assemble_volume_of_fluid_system(volume_of_fluid_field, dim-direction-1, update_from_old);
          }
        solve_volume_of_fluid_system (volume_of_fluid_field);
        // Copy current candidate normals.
        // primarily useful for exact linear translation
        sim.solution.block(volume_of_fluidN_block_idx) = sim.old_solution.block(volume_of_fluidN_block_idx);
        update_volume_of_fluid_normals (volume_of_fluid_field, sim.solution);

        sim.current_linearization_point.block(volume_of_fluid_block_idx) = sim.solution.block(volume_of_fluid_block_idx);
        sim.current_linearization_point.block(volume_of_fluidN_block_idx) = sim.solution.block(volume_of_fluidN_block_idx);
      }
    update_volume_of_fluid_composition(advection_field, volume_of_fluid_field, sim.solution);
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::assemble_volume_of_fluid_system (const VolumeOfFluidField<dim> &field,
                                                                   const unsigned int dir,
                                                                   const bool update_from_old)
  {
    TimerOutput::Scope timer (sim.computing_timer, "Assemble volume of fluid system");

    const unsigned int block0_idx = field_struct_for_field_index(0).volume_fraction.block_index;
    const unsigned int block_idx = field.volume_fraction.block_index;

    if (block0_idx!=block_idx)
      {
        // Allocate the system matrix for the current VoF field by
        // reusing the Trilinos sparsity pattern from the matrix stored for
        // composition 0 (this is the place we allocate the matrix at).
        sim.system_matrix.block(block_idx, block_idx).reinit(sim.system_matrix.block(block0_idx, block0_idx));
      }

    sim.system_matrix.block(block_idx, block_idx) = 0;
    sim.system_rhs = 0;

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    const FiniteElement<dim> &volume_of_fluid_fe = (*field.volume_fraction.fe);

    // we have to assemble the term u.grad phi_i * phi_j, which is
    // of total polynomial degree
    //   stokes_deg - 1
    // (or similar for comp_deg). this suggests using a Gauss
    // quadrature formula of order
    //   stokes_deg/2
    // rounded up. do so. (note that x/2 rounded up
    // equals (x+1)/2 using integer division.)
    const unsigned int vof_quadrature_degree = (this->get_parameters().stokes_velocity_degree+1)/2;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     this->get_dof_handler().begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     this->get_dof_handler().end()),
         std::bind (&Assemblers::VolumeOfFluidAssembler<dim>::
                    local_assemble_volume_of_fluid_system,
                    assembler,
                    field,
                    dir,
                    update_from_old,
                    std::placeholders::_1,
                    std::placeholders::_2,
                    std::placeholders::_3),
         std::bind (&VolumeOfFluidHandler<dim>::
                    copy_local_to_global_volume_of_fluid_system,
                    this,
                    std::placeholders::_1),
         internal::Assembly::Scratch::
         VolumeOfFluidSystem<dim> (this->get_fe(),
                                   volume_of_fluid_fe,
                                   this->get_mapping(),
                                   QGauss<dim>(vof_quadrature_degree),
                                   QGauss<dim-1>(vof_quadrature_degree)),
         internal::Assembly::CopyData::
         VolumeOfFluidSystem<dim> (volume_of_fluid_fe));

    sim.system_matrix.compress(VectorOperation::add);
    sim.system_rhs.compress(VectorOperation::add);
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::copy_local_to_global_volume_of_fluid_system (const internal::Assembly::CopyData::VolumeOfFluidSystem<dim> &data)
  {
    // copy entries into the global matrix. note that these local contributions
    // only correspond to the advection dofs, as assembled above
    sim.current_constraints.distribute_local_to_global (data.local_matrix,
                                                        data.local_rhs,
                                                        data.local_dof_indices,
                                                        sim.system_matrix,
                                                        sim.system_rhs);

    /* In the following, we copy DG contributions element by element. This
     * is allowed since there are no constraints imposed on discontinuous fields.
     */
    for (unsigned int f=0; f<GeometryInfo<dim>::max_children_per_face
         * GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if (data.face_contributions_mask[f])
          {
            for (unsigned int i=0; i<data.neighbor_dof_indices[f].size(); ++i)
              {
                sim.system_rhs(data.neighbor_dof_indices[f][i]) += data.local_face_rhs[f][i];
                for (unsigned int j=0; j< data.neighbor_dof_indices[f].size(); ++j)
                  {
                    sim.system_matrix.add(data.neighbor_dof_indices[f][i],
                                          data.neighbor_dof_indices[f][j],
                                          data.local_face_matrices_ext_ext[f](i,j));
                  }
              }
          }
      }
  }



  template <>
  void VolumeOfFluidHandler<2>::update_volume_of_fluid_normals (const VolumeOfFluidField<2> &field,
                                                                LinearAlgebra::BlockVector &solution)
  {
    const unsigned int dim = 2;
    const unsigned int max_degree = 1;

    LinearAlgebra::BlockVector initial_solution;

    TimerOutput::Scope timer (sim.computing_timer, "Reconstruct VolumeOfFluid interfaces");

    initial_solution.reinit(sim.system_rhs, false);

    // Boundary reference
    typename DoFHandler<dim>::active_cell_iterator endc =
      this->get_dof_handler().end ();

    // Number of cells in the local reconstruction stencil
    const unsigned int n_cells_local_stencil = 9;

    Vector<double> local_volume_of_fluids (n_cells_local_stencil);
    std::vector<Point<dim>> stencil_unit_cell_centers (n_cells_local_stencil);
    std::vector<typename DoFHandler<dim>::active_cell_iterator> neighbor_cells(n_cells_local_stencil);

    const unsigned int n_candidate_normals_per_dim = 3; // Named variable for number of candidate sums for each dimension
    std::vector<double> strip_sums (dim * n_candidate_normals_per_dim);

    // Named variable for number of candidate interface normal vectors for the reconstruction
    const unsigned int n_candidate_normals = dim*n_candidate_normals_per_dim+1;
    std::vector<Tensor<1, dim, double>> normals (n_candidate_normals);
    std::vector<double> d_vals (n_candidate_normals);
    std::vector<double> errs (n_candidate_normals);

    // Variables to do volume calculations

    QGauss<dim> quadrature(max_degree);

    std::vector<double> xFEM_values(quadrature.size());

    const FiniteElement<dim> &system_fe = this->get_fe();

    FEValues<dim> fevalues(this->get_mapping(), system_fe, quadrature,
                           update_JxW_values);

    // Normal holding vars
    Point<dim> reconstruction_stencil_unit_cell_center;
    Tensor<1, dim, double> normal;
    double d;

    for (unsigned int i=0; i<dim; ++i)
      reconstruction_stencil_unit_cell_center[i] = 0.5;

    std::vector<types::global_dof_index> cell_dof_indicies (system_fe.dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices (system_fe.dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int volume_of_fluid_c_index = volume_of_fluid_var.first_component_index;
    const unsigned int volume_of_fluid_ind
      = this->get_fe().component_to_system_index(volume_of_fluid_c_index, 0);

    const FEVariable<dim> &volume_of_fluidN_var = field.reconstruction;
    const unsigned int volume_of_fluidN_c_index = volume_of_fluidN_var.first_component_index;
    const unsigned int volume_of_fluidN_blockidx = volume_of_fluidN_var.block_index;

    const FEVariable<dim> &volume_of_fluidLS_var = field.level_set;
    const unsigned int volume_of_fluidLS_c_index = volume_of_fluidLS_var.first_component_index;
    const unsigned int n_volume_of_fluidLS_dofs = volume_of_fluidLS_var.fe->dofs_per_cell;
    const unsigned int volume_of_fluidLS_blockidx = volume_of_fluidLS_var.block_index;

    //Iterate over cells
    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Obtain data for this cell and neighbors
        cell->get_dof_indices (local_dof_indices);
        const double cell_volume_of_fluid = solution(local_dof_indices[volume_of_fluid_ind]);

        normal[0] = 0.0;
        normal[1] = 0.0;
        d = -1.0;

        if (cell_volume_of_fluid > 1.0 - volume_fraction_threshold)
          {
            d = 1.0;
            initial_solution(local_dof_indices[volume_of_fluid_ind]) = 1.0;
          }
        else if (cell_volume_of_fluid < volume_fraction_threshold)
          {
            d = -1.0;
            initial_solution(local_dof_indices[volume_of_fluid_ind]) = 0.0;
          }
        else
          {
            initial_solution(local_dof_indices[volume_of_fluid_ind]) = cell_volume_of_fluid;

            // Get references to neighboring cells to build stencil references
            //
            // Due to mesh structure, we need to obtain the cells by
            // considering the pattern of neighboring cells
            //
            // Indicies used for the stencil references follow the
            // pattern
            //
            // 6 7 8
            // 3 4 5
            // 0 1 2
            const unsigned int stencil_side_cell_count = 3;
            for (unsigned int i = 0; i < stencil_side_cell_count; ++i)
              {
                typename DoFHandler<dim>::active_cell_iterator cen; // holding variable for cell at column center
                if (i == 0 || i == 2)
                  {
                    // Not on center column, so obtain the center cell of the appropriate row
                    const unsigned int neighbor_no = (i/2);
                    const typename DoFHandler<dim>::face_iterator face = cell->face (neighbor_no);
                    if ((face->at_boundary() && !cell->has_periodic_neighbor(neighbor_no)) ||
                        face->has_children())
                      cen = endc;
                    else
                      {
                        const typename DoFHandler<dim>::cell_iterator neighbor =
                          cell->neighbor_or_periodic_neighbor(neighbor_no);
                        if (neighbor->level() == cell->level() &&
#if DEAL_II_VERSION_GTE(9,2,0)
                            neighbor->is_active())
#else
                            neighbor->active())
#endif
                          cen = neighbor;
                        else
                          cen = endc;
                      }
                  }
                else
                  {
                    // On center column, so current cell is center of column
                    cen = cell;
                  }

                for (unsigned int j = 0; j < stencil_side_cell_count; ++j)
                  {
                    typename DoFHandler<dim>::active_cell_iterator curr; // Variable for cell at current stencil location
                    if (cen == endc)
                      {
                        // Current column center is not in mesh, so assume that
                        // the desired piece of the stencil is also outside.
                        curr = endc;
                      }
                    else
                      {
                        // Current column center exists, so obtain correct cell reference
                        if (j == 0 || j == 2)
                          {
                            //
                            const unsigned int neighbor_no = 2+(j/2);
                            const typename DoFHandler<dim>::face_iterator face = cen->face (neighbor_no);
                            if ((face->at_boundary() && !cen->has_periodic_neighbor(neighbor_no)) ||
                                face->has_children())
                              curr = endc;
                            else
                              {
                                const typename DoFHandler<dim>::cell_iterator neighbor =
                                  cen->neighbor_or_periodic_neighbor(neighbor_no);
                                if (neighbor->level() == cell->level() &&
#if DEAL_II_VERSION_GTE(9,2,0)
                                    neighbor->is_active())
#else
                                    neighbor->active())
#endif
                                  curr = neighbor;
                                else
                                  curr = endc;
                              }
                          }
                        else
                          {
                            // Current stencil reference is column center column center
                            curr = cen;
                          }
                      }
                    if (curr != endc)
                      {
                        // Cell reference is valid, so get data
                        curr->get_dof_indices (cell_dof_indicies);
                        stencil_unit_cell_centers[3 * j + i] = Point<dim> (-1.0 + i,
                                                                           -1.0 + j);
                      }
                    else
                      {
                        // Cell reference is invalid, so replace with current
                        // cell to reduce branching complexity in the later
                        // algorithm
                        cell->get_dof_indices (cell_dof_indicies);
                        stencil_unit_cell_centers[3 * j + i] = Point<dim> (0.0,
                                                                           0.0);
                      }
                    local_volume_of_fluids (3 * j + i) = solution (cell_dof_indicies[volume_of_fluid_ind]);
                    neighbor_cells[3 * j + i] = curr;
                  }
              }

            // Gather cell strip sums
            //
            // Sums are indexed as
            // [n_sums_per_dim * parallel_dim + ind]
            const unsigned int n_sums_per_dim = 3;
            for (unsigned int i = 0; i < dim * n_candidate_normals_per_dim; ++i)
              strip_sums[i] = 0.0;

            for (unsigned int i = 0; i < stencil_side_cell_count; ++i)
              {
                for (unsigned int j = 0; j < stencil_side_cell_count; ++j)
                  {
                    strip_sums[n_sums_per_dim * 0 + i] += local_volume_of_fluids (stencil_side_cell_count * j + i);
                    strip_sums[n_sums_per_dim * 1 + j] += local_volume_of_fluids (stencil_side_cell_count * j + i);
                  }
              }

            // Calculate normal vectors for the 6 candidates from the efficient
            // least squares approach
            //
            // Labeling the sums of the perpendicular strips as
            // 0 1 2
            // L C R
            //
            // For each dimension consider the interface normals implied by the
            // 3 divided differences
            //
            // L-C/1
            // C-R/2
            // L-R/1
            //
            for (unsigned int di = 0; di < dim; ++di)
              {
                // Get index other dimension
                unsigned int di2 = (di == 0) ? 1 : 0;
                for (unsigned int i = 0; i < 3; ++i)
                  {
                    normals[n_candidate_normals_per_dim * di + i][di] = 0.0;
                    normals[n_candidate_normals_per_dim * di + i][di2] = 0.0;
                    if (i==0 || i == 2)
                      {
                        // Positive sum in difference is L
                        normals[n_candidate_normals_per_dim * di + i][di] += strip_sums[n_sums_per_dim * di + 0];
                        normals[3 * di + i][di2] += 1.0;
                      }
                    else
                      {
                        // Positive sum in difference is C
                        normals[3 * di + i][di] += strip_sums[3 * di + 1];
                        normals[3 * di + i][di2] += 0.0;
                      }
                    if (i == 0)
                      {
                        // Negative sum in difference is C
                        normals[3 * di + i][di] -= strip_sums[3 * di + 1];
                        normals[3 * di + i][di2] += 0.0;
                      }
                    else
                      {
                        // Negative sum in difference is R
                        normals[3 * di + i][di] -= strip_sums[3 * di + 2];
                        normals[3 * di + i][di2] += 1.0;
                      }

                    if (strip_sums[3 * di2 + 2] > strip_sums[3 * di2 + 0])
                      {
                        // There is more fluid in in area above the interface on the stencil,
                        // so flip normal direction
                        normals[3 * di + i][di2] *= -1.0;
                      }
                  }
              }

            // Add time extrapolated local normal as candidate
            // this is not expected to be the best candidate in general, but
            // should result in exact reconstruction for linear interface
            // translations
            // Inclusion of this candidate will not reduce accuracy due to it
            // only being selected if it produces a better interface
            // approximation than the ELS candidates. Note that this will
            // render linear translation problems less dependent on the
            // interface reconstruction, so other tests will also be necessary.
            for (unsigned int i=0; i<dim; ++i)
              normals[6][i] = solution(local_dof_indices[system_fe
                                                         .component_to_system_index(volume_of_fluidN_c_index+i, 0)]);

            // If candidate normal too small, remove from consideration
            if (normals[6]*normals[6]< volume_fraction_threshold)
              {
                normals[6][0] = 0;
                normals[6][1] = 0;
              }

            unsigned int index_of_best_normal = 0;
            {
              fevalues.reinit(cell);
              const std::vector<double> weights = fevalues.get_JxW_values();

              double cell_vol = 0.0;
              for (unsigned int j=0; j<weights.size(); ++j)
                {
                  cell_vol+=weights[j];
                }
              for (unsigned int nind = 0; nind < n_candidate_normals; ++nind)
                {
                  errs[nind] = 0.0;
                  const double normal_norm = normals[nind].norm_square();

                  if (normal_norm > volume_fraction_threshold) // If candidate normal too small set error to maximum
                    {
                      d_vals[nind] = VolumeOfFluid::Utilities::compute_interface_location_newton<dim> (
                                       max_degree,
                                       normals[nind],
                                       cell_volume_of_fluid,
                                       cell_vol,
                                       volume_of_fluid_reconstruct_epsilon,
                                       quadrature.get_points(), weights);
                    }
                  else
                    {
                      errs[nind] = 9.0;
                    }
                }
            }

            for (unsigned int i = 0; i < n_cells_local_stencil; ++i)
              {
                if (neighbor_cells[i] == endc)
                  {
                    continue;
                  }

                fevalues.reinit(neighbor_cells[i]);

                const std::vector<double> weights = fevalues.get_JxW_values();

                double cell_vol = 0.0;
                for (unsigned int j=0; j<weights.size(); ++j)
                  {
                    cell_vol+=weights[j];
                  }

                for (unsigned int nind = 0; nind < n_candidate_normals; ++nind)
                  {
                    const double normal_norm = normals[nind]*normals[nind];

                    if (normal_norm > volume_fraction_threshold) // If candidate normal too small skip as set to max already
                      {
                        double dot = 0.0;
                        for (unsigned int di = 0; di < dim; ++di)
                          dot += normals[nind][di] * stencil_unit_cell_centers[i][di];
                        const double n_volume_of_fluid = VolumeOfFluid::Utilities::compute_fluid_volume<dim> (max_degree, normals[nind], d_vals[nind]-dot,
                                                         quadrature.get_points(), weights)/cell_vol;
                        const double cell_err = local_volume_of_fluids (i) - n_volume_of_fluid;
                        errs[nind] += cell_err * cell_err;
                      }
                  }
              }

            for (unsigned int nind = 0; nind < n_candidate_normals; ++nind)
              {
                if (errs[index_of_best_normal] >= errs[nind])
                  index_of_best_normal = nind;
              }

            normal = normals[index_of_best_normal];
            d = d_vals[index_of_best_normal];
          }

        for (unsigned int i=0; i<dim; ++i)
          initial_solution (local_dof_indices[system_fe
                                              .component_to_system_index(volume_of_fluidN_c_index+i, 0)]) = normal[i];

        initial_solution (local_dof_indices[system_fe
                                            .component_to_system_index(volume_of_fluidN_c_index+dim, 0)]) = d;

        for (unsigned int i=0; i<n_volume_of_fluidLS_dofs; ++i)
          {
            // Recenter unit cell on origin
            Tensor<1, dim, double> recentered_support_point = volume_of_fluidLS_var.fe->unit_support_point(i)-reconstruction_stencil_unit_cell_center;
            initial_solution (local_dof_indices[system_fe
                                                .component_to_system_index(volume_of_fluidLS_c_index, i)])
              = d-recentered_support_point*normal;
          }
      }

    initial_solution.compress(VectorOperation::insert);

    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    solution.block(volume_of_fluidN_blockidx) = initial_solution.block(volume_of_fluidN_blockidx);
    solution.block(volume_of_fluidLS_blockidx) = initial_solution.block(volume_of_fluidLS_blockidx);
  }



  template <>
  void VolumeOfFluidHandler<3>::update_volume_of_fluid_normals (const VolumeOfFluidField<3> &/*field*/,
                                                                LinearAlgebra::BlockVector &/*solution*/)
  {
    Assert(false, ExcNotImplemented());
  }



  template <>
  void VolumeOfFluidHandler<2>::update_volume_of_fluid_composition (const typename Simulator<2>::AdvectionField &composition_field,
                                                                    const VolumeOfFluidField<2> &volume_of_fluid_field,
                                                                    LinearAlgebra::BlockVector &solution)
  {
    const unsigned int dim = 2;

    LinearAlgebra::BlockVector initial_solution;

    TimerOutput::Scope timer (sim.computing_timer, "Compute VolumeOfFluid compositions");

    initial_solution.reinit(sim.system_rhs, false);

    // Normal holding vars
    Point<dim> reconstruction_stencil_unit_cell_center;

    for (unsigned int i=0; i<dim; ++i)
      reconstruction_stencil_unit_cell_center[i] = 0.5;

    const FiniteElement<dim> &system_fe = this->get_fe();

    std::vector<types::global_dof_index> local_dof_indices (system_fe.dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = volume_of_fluid_field.volume_fraction;
    const unsigned int volume_of_fluid_c_index = volume_of_fluid_var.first_component_index;
    const unsigned int volume_of_fluid_ind
      = system_fe.component_to_system_index(volume_of_fluid_c_index, 0);

    const FEVariable<dim> &volume_of_fluidN_var = volume_of_fluid_field.reconstruction;
    const unsigned int volume_of_fluidN_c_index = volume_of_fluidN_var.first_component_index;

    const unsigned int base_element = composition_field.base_element(this->introspection());
    const std::vector<Point<dim> > support_points = system_fe.base_element(base_element).get_unit_support_points();

    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned())
          continue;

        cell->get_dof_indices (local_dof_indices);
        const double cell_volume_of_fluid = solution(local_dof_indices[volume_of_fluid_ind]);

        Tensor<1, dim, double> normal;

        for (unsigned int i=0; i<dim; ++i)
          normal[i] = solution(local_dof_indices[system_fe
                                                 .component_to_system_index(volume_of_fluidN_c_index+i, 0)]);

        // Compute L1 normalized normal vector
        // (Use of L1 normalization makes calculating the correct correction
        // factor much simpler to retain $0\leq C\leq 1$ bound much simpler)
        double normal_l1_norm = 0.0;
        for (unsigned int i=0; i<dim; ++i)
          {
            normal_l1_norm += std::abs(normal[i]);
          }
        //Calculate correct factor to retain vol frac and [0,1] bound
        double correction_factor = (normal_l1_norm<volume_fraction_threshold)
                                   ?
                                   0.0
                                   :
                                   2.0*(0.5-abs(cell_volume_of_fluid-0.5))/normal_l1_norm;
        for (unsigned int i=0; i<system_fe.base_element(base_element).dofs_per_cell; ++i)
          {
            const unsigned int system_local_dof
              = system_fe.component_to_system_index(composition_field.component_index(sim.introspection),
                                                    /*dof index within component*/i);

            Tensor<1, dim, double> recentered_support_point = support_points[i]-reconstruction_stencil_unit_cell_center;

            const double value = cell_volume_of_fluid - correction_factor*(recentered_support_point*normal);

            initial_solution(local_dof_indices[system_local_dof]) = value;
          }

      }


    const unsigned int blockidx = composition_field.block_index(this->introspection());
    solution.block(blockidx) = initial_solution.block(blockidx);
  }



  template <>
  void VolumeOfFluidHandler<3>::update_volume_of_fluid_composition (const typename Simulator<3>::AdvectionField &/*composition_field*/,
                                                                    const VolumeOfFluidField<3> &/*volume_of_fluid_field*/,
                                                                    LinearAlgebra::BlockVector &/*solution*/)
  {
    Assert(false, ExcNotImplemented());
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::set_initial_volume_fractions ()
  {
    for (unsigned int f=0; f<n_volume_of_fluid_fields; ++f)
      {

        switch (initialization_data_type[f])
          {
            case VolumeOfFluid::VolumeOfFluidInputType::composition:
              initialize_from_composition_field (data[f]);
              break;
            case VolumeOfFluid::VolumeOfFluidInputType::level_set:
              initialize_from_level_set (data[f]);
              break;
            default:
              Assert(false, ExcNotImplemented ());
          }

        const unsigned int volume_of_fluidN_blockidx = data[f].reconstruction.block_index;
        const unsigned int volume_of_fluidLS_blockidx = data[f].level_set.block_index;
        update_volume_of_fluid_normals (data[f], sim.solution);
        sim.old_solution.block(volume_of_fluidN_blockidx) = sim.solution.block(volume_of_fluidN_blockidx);
        sim.old_old_solution.block(volume_of_fluidN_blockidx) = sim.solution.block(volume_of_fluidN_blockidx);
        sim.old_solution.block(volume_of_fluidLS_blockidx) = sim.solution.block(volume_of_fluidLS_blockidx);
        sim.old_old_solution.block(volume_of_fluidLS_blockidx) = sim.solution.block(volume_of_fluidLS_blockidx);

        // Update associated composition field
        const typename Simulator<dim>::AdvectionField composition_field = Simulator<dim>::AdvectionField::composition(data[f].composition_index);
        update_volume_of_fluid_composition (composition_field, data[f], sim.solution);
        const unsigned int volume_of_fluid_C_blockidx = composition_field.block_index(this->introspection());
        sim.old_solution.block(volume_of_fluid_C_blockidx) = sim.solution.block(volume_of_fluid_C_blockidx);
        sim.old_old_solution.block(volume_of_fluid_C_blockidx) = sim.solution.block(volume_of_fluid_C_blockidx);
      }
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::initialize_from_composition_field (
    const VolumeOfFluidField<dim> &field)
  {
    LinearAlgebra::BlockVector initial_solution;

    initial_solution.reinit(sim.system_rhs, false);

    const QIterated<dim> quadrature (QMidpoint<1>(), n_init_samples);
    FEValues<dim, dim> fe_init (this->get_mapping(), this->get_fe(), quadrature,
                                update_JxW_values | update_quadrature_points);

    std::vector<types::global_dof_index>
    local_dof_indices (this->get_fe().dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int component_index = volume_of_fluid_var.first_component_index;
    const unsigned int blockidx = volume_of_fluid_var.block_index;
    const unsigned int volume_of_fluid_ind
      = this->get_fe().component_to_system_index(component_index, 0);

    // Initialize state based on provided function
    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Calculate approximation for volume
        cell->get_dof_indices (local_dof_indices);

        fe_init.reinit (cell);

        double volume_of_fluid_val = 0.0;
        double cell_vol = 0.0;

        for (unsigned int i = 0; i < fe_init.n_quadrature_points; ++i)
          {
            const double fraction_at_point = this->get_initial_composition_manager().initial_composition(fe_init.quadrature_point(i),
                                             field.composition_index);
            volume_of_fluid_val += fraction_at_point * fe_init.JxW (i);
            cell_vol += fe_init.JxW(i);
          }

        volume_of_fluid_val /= cell_vol;

        volume_of_fluid_val = std::min(volume_of_fluid_val, 1.0);
        volume_of_fluid_val = std::max(volume_of_fluid_val, 0.0);

        initial_solution (local_dof_indices[volume_of_fluid_ind]) = volume_of_fluid_val;
      }

    initial_solution.compress(VectorOperation::insert);

    // Apply constraints and update solution blocks. This is duplicated from
    // Simulator<dim>::set_initial_temperature_and_compositional_fields()
    // in order to separate the volume of fluid algorithm as much as possible from
    // the rest of the code.
    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    sim.solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_old_solution.block(blockidx) = initial_solution.block(blockidx);
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::initialize_from_level_set (
    const VolumeOfFluidField<dim> &field)
  {
    LinearAlgebra::BlockVector initial_solution;

    initial_solution.reinit(sim.system_rhs, false);

    const QIterated<dim> quadrature (QMidpoint<1>(), n_init_samples);
    FEValues<dim, dim> fe_init (this->get_mapping(),
                                this->get_fe(),
                                quadrature,
                                update_JxW_values | update_quadrature_points);

    const double h = 1.0/n_init_samples;

    std::vector<types::global_dof_index>
    local_dof_indices (this->get_fe().dofs_per_cell);

    const FEVariable<dim> &volume_of_fluid_var = field.volume_fraction;
    const unsigned int component_index = volume_of_fluid_var.first_component_index;
    const unsigned int blockidx = volume_of_fluid_var.block_index;
    const unsigned int volume_of_fluid_ind
      = this->get_fe().component_to_system_index(component_index, 0);

    // Initialize state based on provided function
    for (auto cell : this->get_dof_handler().active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Calculate approximation for volume
        cell->get_dof_indices (local_dof_indices);

        const double cell_diam = cell->diameter();
        const double d_func = this->get_initial_composition_manager().initial_composition(cell->barycenter(),
                              field.composition_index);
        fe_init.reinit (cell);

        double volume_of_fluid_val = 0.0;
        double cell_vol = 0.0;

        if (d_func <=-0.5*cell_diam)
          {
            volume_of_fluid_val = 0.0;
          }
        else if (d_func >= 0.5*cell_diam)
          {
            volume_of_fluid_val = 1.0;
          }
        else
          {

            // For each quadrature point compute an approximation to the fluid fraction in the surrounding region
            for (unsigned int i = 0; i < fe_init.n_quadrature_points; ++i)
              {
                double d = 0.0;
                Tensor<1, dim, double> grad;
                Point<dim> xU = quadrature.point (i);

                // Get an approximation to local normal at the closest interface (level set gradient)
                // and the distance to the closest interface (value of level set function)
                for (unsigned int di = 0; di < dim; ++di)
                  {
                    Point<dim> xH, xL;
                    xH = xU;
                    xL = xU;
                    xH[di] += 0.5*h;
                    xL[di] -= 0.5*h;
                    const double dH = this->get_initial_composition_manager().initial_composition(cell->intermediate_point(xH),
                                                                                                  field.composition_index);
                    const double dL = this->get_initial_composition_manager().initial_composition(cell->intermediate_point(xL),
                                                                                                  field.composition_index);
                    grad[di] = (dL-dH);
                    d += (0.5/dim)*(dH+dL);
                  }
                // Use the basic fluid fraction formula to compute an approximation to the fluid fraction
                const double fraction_at_point = VolumeOfFluid::Utilities::compute_fluid_fraction (grad, d);
                volume_of_fluid_val += fraction_at_point * fe_init.JxW (i);
                cell_vol += fe_init.JxW (i);
              }
            volume_of_fluid_val /= cell_vol;
          }

        initial_solution (local_dof_indices[volume_of_fluid_ind]) = volume_of_fluid_val;
      }

    initial_solution.compress(VectorOperation::insert);

    // Apply constraints and update solution blocks. This is duplicated from
    // Simulator<dim>::set_initial_temperature_and_compositional_fields()
    // in order to separate the volume of fluid algorithm as much as possible from
    // the rest of the code.
    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    sim.solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_solution.block(blockidx) = initial_solution.block(blockidx);
    sim.old_old_solution.block(blockidx) = initial_solution.block(blockidx);
  }



  template <int dim>
  void VolumeOfFluidHandler<dim>::solve_volume_of_fluid_system (const VolumeOfFluidField<dim> &field)
  {
    const unsigned int block_idx = field.volume_fraction.block_index;

    TimerOutput::Scope timer (sim.computing_timer, "Solve volume of fluid system");
    this->get_pcout() << "   Solving volume of fluid system... " << std::flush;

    const double tolerance = std::max(1e-50,
                                      volume_of_fluid_solver_tolerance*sim.system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);

#ifdef ASPECT_USE_PETSC
    SolverCG<LinearAlgebra::Vector> solver(solver_control);
    LinearAlgebra::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));
#else
    TrilinosWrappers::SolverCG solver(solver_control);
    TrilinosWrappers::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));
#endif

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because only have a ConstraintMatrix
    // for the whole system, current_linearization_point contains our initial guess.
    LinearAlgebra::BlockVector distributed_solution (
      this->introspection().index_sets.system_partitioning,
      this->get_mpi_communicator());
    distributed_solution.block(block_idx) = sim.current_linearization_point.block (block_idx);

    sim.current_constraints.set_zero(distributed_solution);

    // solve the linear system:
    try
      {
        solver.solve (sim.system_matrix.block(block_idx,block_idx),
                      distributed_solution.block(block_idx),
                      sim.system_rhs.block(block_idx),
                      precondition);
      }
    // if the solver fails, report the error from processor 0 with some additional
    // information about its location, and throw a quiet exception on all other
    // processors
    catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          AssertThrow (false,
                       ExcMessage (std::string("The iterative advection solver "
                                               "did not converge. It reported the following error:\n\n")
                                   +
                                   exc.what()))
          else
            throw QuietException();
      }

    sim.current_constraints.distribute (distributed_solution);
    sim.solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    this->get_pcout() << solver_control.last_step()
                      << " iterations." << std::endl;

    // Do not add VolumeOfFluid solver iterations to statistics, duplication due to
    // dimensional splitting results in incorrect line formatting (lines of
    // data split inconsistently with missing values)
  }

}



namespace aspect
{
#define INSTANTIATE(dim) \
  template class VolumeOfFluidHandler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
