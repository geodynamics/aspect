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

    this->get_signals().edit_finite_element_variables.connect(
      [&](std::vector<VariableDeclaration<dim> > &vars)
    {
      this->edit_finite_element_variables(vars);
    });

    this->get_signals().post_set_initial_state.connect(
      [&](const SimulatorAccess<dim> &)
    {
      this->set_initial_volume_fractions();
    });
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
      prm.declare_entry("Volume of fluid initialization type", "",
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
      const std::vector<std::string> x_initialization_type = Utilities::split_string_list(prm.get("Volume of fluid initialization type"));

      initialization_data_type = std::vector<VolumeOfFluid::VolumeOfFluidInputType::Kind> (n_volume_of_fluid_fields,
                                 VolumeOfFluid::VolumeOfFluidInputType::composition);

      for (const auto &p : x_initialization_type)
        {
          // each entry has the format (white space is optional):
          // <name> : <value (might have spaces)>
          //
          // first tease apart the two halves
          const std::vector<std::string> split_parts = Utilities::split_string_list (p, ':');
          AssertThrow (split_parts.size() == 2,
                       ExcMessage ("The format for "
                                   "<Initial composition model/Volume of fluid initialization method> "
                                   "volume of fluid initialization met "
                                   "requires that each entry " "has the form "
                                   "`<name of field> : <method>', but this "
                                   "does not match the number of colons in the "
                                   "entry <" + p + ">."));

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
                                   "<Initial composition model/Volume of fluid initialization method>, but "
                                   "there is no field with this name."));

          const unsigned int compositional_field_index = std::distance(names_of_compositional_fields.begin(),
                                                                       field_name_iterator);

          AssertThrow (compositional_field_methods[compositional_field_index]
                       == Parameters<dim>::AdvectionFieldMethod::volume_of_fluid,
                       ExcMessage ("The field <" + key + "> appears in the parameter "
                                   "<Initial composition model/Volume of fluid initialization method>, "
                                   "but is not advected by a particle method."));

          AssertThrow (std::count(names_of_compositional_fields.begin(),
                                  names_of_compositional_fields.end(), key) == 1,
                       ExcMessage ("Name of field <" + key + "> appears more "
                                   "than once in the parameter "
                                   "<Initial composition model/Volume of fluid initialization type>."));

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
                ExcMessage("Volume of Fluid Interface Tracking currently assumes incompressibility."));

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

    using CellFilter = FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>;

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
}



namespace aspect
{
#define INSTANTIATE(dim) \
  template class VolumeOfFluidHandler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
