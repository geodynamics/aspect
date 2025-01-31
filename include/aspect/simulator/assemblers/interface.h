/*
  Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_simulator_assemblers_interface_h
#define _aspect_simulator_assemblers_interface_h

#include <aspect/global.h>
#include <aspect/heating_model/interface.h>
#include <aspect/material_model/interface.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  template <int dim>
  class Simulator;

  struct AdvectionField;

  /**
   * A namespace that is used for internal scratch objects, i.e. objects that define
   * the inputs and outputs for the individual assembler objects. The different classes
   * are provided for the different systems of equations, and the inputs are filled
   * by the <code>local_assemble_...</code> functions in <code>assembly.cc</code>,
   * before being handed over to the assemblers.
   */
  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        /**
         * Scratch objects are used to store information about a cell that is
         * necessary for assembling matrix and right hand side terms for this
         * cell. The ScratchBase class acts as a empty base class for
         * individual scratch objects for the different equations.
         */
        template <int dim>
        struct ScratchBase
        {
          ScratchBase()
            :
            cell(),
            face_number(numbers::invalid_unsigned_int)
          {}

          ScratchBase(const ScratchBase &scratch)
            :
            cell(scratch.cell),
            face_number(scratch.face_number)
          {}

          virtual ~ScratchBase ()  = default;

          /**
           * Cell object on which we currently operate.
           */
          typename DoFHandler<dim>::active_cell_iterator cell;

          /**
           * The number of the face object with respect to the current
           * cell on which we operate. If we currently
           * operate on a cell, this member is set to
           * numbers::invalid_unsigned_int.
           */
          unsigned face_number;
        };

        /**
         * A scratch object to store all necessary information to assemble
         * the Stokes preconditioner terms.
         */
        template <int dim>
        struct StokesPreconditioner: public ScratchBase<dim>
        {
          StokesPreconditioner (const FiniteElement<dim> &finite_element,
                                const Quadrature<dim>    &quadrature,
                                const Mapping<dim>       &mapping,
                                const UpdateFlags         update_flags,
                                const unsigned int        n_compositional_fields,
                                const unsigned int        stokes_dofs_per_cell,
                                const bool                add_compaction_pressure,
                                const bool                rebuild_matrix,
                                const bool                use_bfbt);
          StokesPreconditioner (const StokesPreconditioner &scratch);

          ~StokesPreconditioner () override;

          FEValues<dim> finite_element_values;

          void reinit (const typename DoFHandler<dim>::active_cell_iterator &cell_ref);

          std::vector<types::global_dof_index> local_dof_indices;
          std::vector<unsigned int>            dof_component_indices;
          std::vector<SymmetricTensor<2,dim>> grads_phi_u;
          std::vector<double>                  div_phi_u;
          std::vector<double>                  phi_p;
          std::vector<Tensor<1,dim>>           phi_u;
          std::vector<double>                  phi_p_c;
          std::vector<Tensor<1,dim>>          grad_phi_p;

          /**
           * Material model inputs and outputs computed at the current
           * linearization point.
           */
          MaterialModel::MaterialModelInputs<dim> material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> material_model_outputs;

          /**
           * Whether the Stokes matrix should be rebuild during this
           * assembly. If the matrix does not change, assembling the right
           * hand side is sufficient.
           */
          const bool rebuild_stokes_matrix;
        };



        /**
         * A scratch object to store all necessary information to assemble
         * the terms in the Stokes equations.
         * We derive the StokesSystem scratch class from the
         * StokesPreconditioner class, because all the objects that
         * are necessary for the assembly of the preconditioner are also
         * needed for the actual system matrix and right hand side, plus some
         * extra data that we need for the time stepping and traction boundaries
         * on the right hand side.
         */
        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &finite_element,
                        const Mapping<dim>       &mapping,
                        const Quadrature<dim>    &quadrature,
                        const Quadrature<dim-1>  &face_quadrature,
                        const UpdateFlags         update_flags,
                        const UpdateFlags         face_update_flags,
                        const unsigned int        n_compositional_fields,
                        const unsigned int        stokes_dofs_per_cell,
                        const bool                add_compaction_pressure,
                        const bool                use_reference_density_profile,
                        const bool                rebuild_stokes_matrix,
                        const bool                rebuild_newton_stokes_matrix,
                        const bool                use_bfbt);

          StokesSystem (const StokesSystem<dim> &scratch);

          FEFaceValues<dim> face_finite_element_values;

          using StokesPreconditioner<dim>::reinit;

          void reinit (const typename DoFHandler<dim>::active_cell_iterator &cell_ref,
                       const unsigned face_number_ref);

          std::vector<Tensor<1,dim>>          phi_u;
          std::vector<Tensor<1,dim>>          velocity_values;
          std::vector<double>                  velocity_divergence;
          std::vector<Tensor<1,dim>>          temperature_gradients;

          /**
           * Material model inputs and outputs computed at the current
           * linearization point.
           *
           * In contrast to the variables above, the following two
           * variables are used in the assembly at quadrature points
           * on faces, not on cells.
           */
          MaterialModel::MaterialModelInputs<dim> face_material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> face_material_model_outputs;

          /**
           * In some approximations of the Stokes equations the density used
           * for the mass conservation (/continuity) equation is some form
           * of reference density, while the density used for calculating
           * the buoyancy force is the full density. In case such a formulation
           * is used the reference density (and its derivative in depth
           * direction) is queried from the adiabatic conditions plugin
           * and is stored in these variables.
           */
          std::vector<double> reference_densities;
          std::vector<double> reference_densities_depth_derivative;

          /**
           * Whether the Newton solver Stokes matrix should be rebuild during
           * this assembly. If the matrix does not change, assembling the right
           * hand side is sufficient.
           */
          const bool rebuild_newton_stokes_matrix;
        };



        /**
         * A scratch object to store all necessary information to assemble
         * the terms in the advection equations.
         */
        template <int dim>
        struct AdvectionSystem: public ScratchBase<dim>
        {
          AdvectionSystem (const FiniteElement<dim> &finite_element,
                           const FiniteElement<dim> &advection_element,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &quadrature,
                           const Quadrature<dim-1>  &face_quadrature,
                           const UpdateFlags         update_flags,
                           const UpdateFlags         face_update_flags,
                           const unsigned int        n_compositional_fields,
                           const typename Simulator<dim>::AdvectionField     &field);
          AdvectionSystem (const AdvectionSystem &scratch);

          FEValues<dim> finite_element_values;

          void reinit (const typename DoFHandler<dim>::active_cell_iterator &cell_ref);

          std::unique_ptr<FEFaceValues<dim>>    face_finite_element_values;
          std::unique_ptr<FEFaceValues<dim>>    neighbor_face_finite_element_values;
          std::unique_ptr<FESubfaceValues<dim>> subface_finite_element_values;

          std::vector<types::global_dof_index>   local_dof_indices;

          /**
           * Variables describing the values and gradients of the
           * shape functions at the quadrature points, as they are
           * used in the advection assembly function. note that the sizes
           * of these arrays are equal to the number of shape functions
           * corresponding to the currently advected field (and not all of the
           * existing fields), and that they are also correspondingly indexed.
           */
          std::vector<double>         phi_field;
          std::vector<Tensor<1,dim>> grad_phi_field;
          std::vector<double>         laplacian_phi_field;
          std::vector<double>         face_phi_field;
          std::vector<Tensor<1,dim>> face_grad_phi_field;
          std::vector<double>         neighbor_face_phi_field;
          std::vector<Tensor<1,dim>> neighbor_face_grad_phi_field;

          std::vector<Tensor<1,dim>> old_velocity_values;
          std::vector<Tensor<1,dim>> old_old_velocity_values;

          std::vector<double>         old_pressure;
          std::vector<double>         old_old_pressure;
          std::vector<Tensor<1,dim>> old_pressure_gradients;
          std::vector<Tensor<1,dim>> old_old_pressure_gradients;

          std::vector<SymmetricTensor<2,dim>> old_strain_rates;
          std::vector<SymmetricTensor<2,dim>> old_old_strain_rates;

          std::vector<double>         old_temperature_values;
          std::vector<double>         old_old_temperature_values;

          std::vector<double>         old_field_values;
          std::vector<double>         old_old_field_values;
          std::vector<Tensor<1,dim>> old_field_grads;
          std::vector<Tensor<1,dim>> old_old_field_grads;
          std::vector<double>         old_field_laplacians;
          std::vector<double>         old_old_field_laplacians;

          std::vector<std::vector<double>> old_composition_values;
          std::vector<std::vector<double>> old_old_composition_values;

          std::vector<double>         current_temperature_values;
          std::vector<Tensor<1,dim>> current_velocity_values;
          std::vector<Tensor<1,dim>> face_current_velocity_values;
          std::vector<Tensor<1,dim>> mesh_velocity_values;
          std::vector<Tensor<1,dim>> face_mesh_velocity_values;

          std::vector<SymmetricTensor<2,dim>> current_strain_rates;
          std::vector<std::vector<double>> current_composition_values;
          std::vector<double>         current_velocity_divergences;

          /**
           * Material model inputs and outputs computed at the current
           * linearization point.
           */
          MaterialModel::MaterialModelInputs<dim> material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> material_model_outputs;

          MaterialModel::MaterialModelInputs<dim> face_material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> face_material_model_outputs;

          MaterialModel::MaterialModelInputs<dim> neighbor_face_material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> neighbor_face_material_model_outputs;

          /**
           * Heating model outputs computed at the quadrature points of the
           * current cell at the time of the current linearization point.
           * As explained in the class documentation of
           * HeatingModel::HeatingModelOutputs each term contains the sum of all
           * enabled heating mechanism contributions.
           */
          HeatingModel::HeatingModelOutputs heating_model_outputs;
          HeatingModel::HeatingModelOutputs face_heating_model_outputs;
          HeatingModel::HeatingModelOutputs neighbor_face_heating_model_outputs;

          /**
           * This pointer contains a struct that can be used to identify the
           * advection field that is currently assembled. It can be used to
           * determine between temperature and the available compositional
           * fields. See the documentation of the AdvectionField class for
           * more details.
           */
          const typename Simulator<dim>::AdvectionField *advection_field;

          /**
           * The amount of entropy viscosity that should be applied to the
           * current cell to stabilize the solution of the advection system.
           */
          double artificial_viscosity;
        };
      }



      /**
       * The CopyData arrays are similar to the Scratch arrays except they are
       * meant as containers for the output of assembler objects. They provide a
       * constructor and some data objects for local matrix, local vectors and
       * the relation between local and global degrees of freedom (a.k.a.
       * <code>local_dof_indices</code>). After all assemblers are finished
       * the objects contain the local contributions of a particular cell to
       * the global matrix and right hand side. This copy data object is then
       * handed over to one of the <code>copy_local_to_global...</code>
       * functions in assembly.cc that copy their content to the global matrix
       * and right hand side vector.
       */
      namespace CopyData
      {
        /**
         * The base class is empty and only allows us to hand over pointers
         * or references of a generic type and later cast them to their actual
         * derived class.
         */
        template <int dim>
        struct CopyDataBase
        {
          virtual ~CopyDataBase () = default;
        };

        /**
         * The Stokes preconditioner object only requires the bare minimum of
         * copy data objects. Matrix contributions and degrees of freedom this
         * cell corresponds to.
         */
        template <int dim>
        struct StokesPreconditioner: public CopyDataBase<dim>
        {
          StokesPreconditioner (const unsigned int stokes_dofs_per_cell);

          StokesPreconditioner (const StokesPreconditioner &data);

          ~StokesPreconditioner () override = default;
          StokesPreconditioner<dim> &operator= (const StokesPreconditioner<dim> &data) = default;

          FullMatrix<double> local_matrix;
          Vector<double> local_inverse_lumped_mass_matrix;
          std::vector<types::global_dof_index> local_dof_indices;

          /**
           * Extract the values listed in @p all_dof_indices only if
           * it corresponds to the Stokes component and copy it to the variable
           * local_dof_indices declared above in the same class as this function
           */
          void extract_stokes_dof_indices(const std::vector<types::global_dof_index> &all_dof_indices,
                                          const Introspection<dim>                   &introspection,
                                          const FiniteElement<dim>                   &finite_element);
        };

        /**
         * Similar to the scratch object the Stokes system requires all
         * data from the Stokes preconditioner copy data class, plus some
         * extras like the right hand side contribution.
         */
        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const unsigned int        stokes_dofs_per_cell,
                        const bool                do_pressure_rhs_compatibility_modification);
          StokesSystem (const StokesSystem<dim> &data);

          ~StokesSystem () override = default;
          StokesSystem<dim> &operator= (const StokesSystem<dim> &data) = default;

          Vector<double> local_rhs;
          Vector<double> local_pressure_shape_function_integrals;
        };

        /**
         * Additionally to the Stokes system the Advection system copy data
         * object also needs to keep track of contributions across faces
         * (mostly for discontinuous elements that contain DG terms).
         */
        template <int dim>
        struct AdvectionSystem: public CopyDataBase<dim>
        {
          /**
           * Constructor.
           *
           * @param finite_element The element that describes the field for
           *    which we are trying to assemble a linear system. <b>Not</b>
           *    the global finite element.
           * @param field_is_discontinuous If true, the field is a DG element.
           */
          AdvectionSystem (const FiniteElement<dim> &finite_element,
                           const bool                field_is_discontinuous);

          /**
           * Local contributions to the global matrix
           * that correspond only to the variables listed in local_dof_indices
           */
          FullMatrix<double>          local_matrix;

          /**
           * Local contributions to the global matrix from the face terms in the
           * discontinuous Galerkin method. These arrays are of a length sufficient
           * to hold one matrix for each possible face or subface of the cell.
           * The discontinuous Galerkin bilinear form contains terms arising from
           * internal (to the cell) values and external (to the cell) values.
           * `_int_ext` and `_ext_int` hold the terms arising from the pairing
           * between a cell and its neighbor, while `_ext_ext` is the pairing
           * of the neighbor's dofs with themselves. In the continuous
           * Galerkin case, these are unused, and set to size zero.
           */
          std::vector<FullMatrix<double>>         local_matrices_int_ext;
          std::vector<FullMatrix<double>>         local_matrices_ext_int;
          std::vector<FullMatrix<double>>         local_matrices_ext_ext;

          /**
           * Local contributions to the right hand side
           * that correspond only to the variables listed in local_dof_indices
           */
          Vector<double>              local_rhs;

          /**
           * Denotes which face matrices have actually been assembled in the DG field
           * assembly. Entries for matrices not used (for example, those corresponding
           * to non-existent subfaces; or faces being assembled by the neighboring cell)
           * are set to false.
           */
          std::vector<bool>               assembled_matrices;

          /**
           * Indices of those degrees of freedom that actually correspond
           * to the temperature or compositional field. since this structure
           * is used to represent just contributions to the advection
           * systems, there will be no contributions to other parts of the
           * system and consequently, we do not need to list here indices
           * that correspond to velocity or pressure degrees (or, in fact
           * any other variable outside the block we are currently considering)
           */
          std::vector<types::global_dof_index>   local_dof_indices;

          /**
           * Indices of the degrees of freedom corresponding to the temperature
           * or composition field on all possible neighboring cells. This is used
           * in the discontinuous Galerkin method. The outer array has a
           * length sufficient to hold one element for each possible face
           * and sub-face of the current cell. The object is not used
           * and has size zero if in the continuous Galerkin case.
           */
          std::vector<std::vector<types::global_dof_index>>   neighbor_dof_indices;
        };
      }
    }
  }

  /**
   * A namespace for the definition of assemblers for the various terms in the
   * linear systems ASPECT solves.
   */
  namespace Assemblers
  {

    /**
     * For a reference cell (which is typically obtained by asking the finite
     * element to be used), determine how many interface matrices are needed.
     * Since interface matrices are needed for as many neighbors as each
     * cell can have, this is the number of faces for the given reference cell
     * times the number of children each of these faces can have. This
     * accommodates the fact that the neighbors of a cell can all be refined,
     * though they can only be refined once.
     */
    unsigned int
    n_interface_matrices (const ReferenceCell &reference_cell);

    /**
     * For a given reference cell, and a given face we are currently
     * assembling on, return which element of an array of size
     * `n_interface_matrices(reference_cell)` to use.
     */
    unsigned int
    nth_interface_matrix (const ReferenceCell &reference_cell,
                          const unsigned int face);

    /**
     * For a given reference cell, and a given face and sub-face we are
     * currently assembling on, return which element of an array of size
     * `n_interface_matrices(reference_cell)` to use.
     */
    unsigned int
    nth_interface_matrix (const ReferenceCell &reference_cell,
                          const unsigned int face,
                          const unsigned int sub_face);

    /**
     * A base class for objects that implement assembly
     * operations.
     *
     * The point of this class is primarily so that we can store
     * pointers to such objects in a list. The objects are created
     * in Simulator::set_assemblers() and destroyed in the destructor of
     * the Simulator object. Derived classes of this base class usually
     * handle groups of terms in the equations that are
     * logically connected (such as all terms in the Stokes equations that
     * appear independent on the selected compressibility formulation). This
     * way selecting a certain set of assembler objects effectively controls
     * which equation is solved.
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Execute this assembler object. This function performs the primary work
         * of an assembler. More precisely, it uses information for the current
         * cell that is stored in @p scratch (like the material properties on
         * this cell and the position of quadrature points) and computes the
         * matrix and right hand side contributions for a set of terms for
         * the given cell. These contributions are stored in @p data. Note, that
         * the data in @p scratch and @p data is shared between all active
         * assemblers so that each assembler should only add contributions to
         * @p data, not overwrite entries in the matrix. After all assemblers
         * have finished, the final content of @p data is distributed into the
         * global matrix and right hand side vector.
         */
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim> &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const = 0;

        /**
         * This function gets called if a MaterialModelOutputs is created
         * and allows the assembler to attach AdditionalOutputs. The
         * function might be called more than once for a
         * MaterialModelOutput, so it is recommended to check if
         * get_additional_output() returns an instance before adding a new
         * one to the additional_outputs vector. By default this function does
         * not create additional outputs.
         *
         * Material models, through functions derived from
         * MaterialModel::Interface::evaluate(), put their computed material
         * parameters into a structure of type MaterialModel::MaterialModelOutputs.
         * By default, material models will compute those parameters that
         * correspond to the member variables of that structure. However,
         * there are situations where parts of the simulator need additional
         * pieces of information; a typical example would be the use of a
         * Newton scheme that also requires the computation of <i>derivatives</i>
         * of material parameters with respect to pressure, temperature, and
         * possibly other variables.
         *
         * The computation of such additional information is controlled by
         * the presence of a collection of pointers in
         * MaterialModel::MaterialModelOutputs that point to additional
         * objects. Whether or not one needs these additional objects depends
         * on what assemblers are selected, or what postprocessing one
         * wants to compute. For the purpose of assembly, the current
         * function creates the additional objects (such as the one that stores
         * derivatives) and adds pointers to them to the collection, based on
         * what this assembler class requires. This function is always called
         * before the material model is evaluated and execute() is called.
         * This ensures the additional material model output is available when
         * execute() is called.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &) const;

        /**
         * A required function for objects that implement the assembly of terms
         * in an equation that requires the computation of residuals
         * (in particular the advection equation in ASPECT).
         * Just like the assemblers itself, the residual
         * that we use to compute the necessary entropy viscosity depend on the
         * equation (i.e. which terms are actually included in the
         * equation). Thus different objects compute different residuals (i.e.
         * the residual for a melt advection equation looks different from the
         * residual for a passive compositional field).
         * For assemblers for the Stokes system, an implementation of this
         * function is not necessary.
         */
        virtual
        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &) const;
    };



    /**
     * A base class for objects that implement assembly
     * operations for advection-diffusion problems.
     *
     * This class implements functions that provide information
     * for stabilization mechanisms.
     */
    template <int dim>
    class AdvectionStabilizationInterface
    {
      public:
        virtual ~AdvectionStabilizationInterface ();

        /**
         * This function returns a representative prefactor for the advection
         * term of the equation for each quadrature point of the current cell.
         * In the non-dimensional case this is simply 1.0, but for other
         * quantities like temperature it is computed using physical units
         * (like density and specific heat capacity).
         * This information is useful for algorithms that depend on the
         * magnitude of individual terms, like stabilization methods.
         */
        virtual
        std::vector<double>
        advection_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const;

        /**
         * This function returns a representative conductivity for the
         * diffusion part of the equation for each quadrature point of the
         * current cell. For the pure advection case this factor is 0.0, but
         * for other quantities like temperature it is
         * computed using physical units (like thermal conductivity).  This
         * information is useful for algorithms depending on the magnitude of
         * individual terms, like stabilization methods.
         */
        virtual
        std::vector<double>
        diffusion_prefactors(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const;
    };



    /**
     * A class that owns member variables representing
     * all assemblers that need to be called when
     * assembling right hand side vectors, matrices, or complete linear
     * systems. We use this approach in order to support the following
     * cases:
     * - Assembling different formulations: When assembling either the
     *   full equations or only the Boussinesq approximation (to give just
     *   two examples), one needs different terms. This could be achieved
     *   using a large number of <code>switch</code> or <code>if</code>
     *   statements in the code, or one could encapsulate each equation
     *   or approximation into a collection of assemblers for this particular
     *   purpose. The approach chosen here in essence allows the
     *   implementation of each set of equations in its own scope, and we
     *   then just need to store a pointer to the object that
     *   assembles the Stokes system (for example) for the selected
     *   approximation. The pointer to this function is stored in the
     *   appropriate member variable of this class.
     * - Sometimes, we want to assemble a number of terms that build on
     *   each other. An example is the addition of free boundary terms
     *   to the Stokes matrix. Rather than having to "know" in one
     *   place about all of the terms that need to be assembled,
     *   we simply add the function that computes these terms as
     *   another object to the appropriate set of assemblers declared
     *   in this class.
     */
    template <int dim>
    class Manager
    {
      public:

        /**
         * Reset the state of the manager and remove all Assemblers.
         */
        void reset ();

        /**
         * A vector of pointers containing all assemblers for the Stokes preconditioner.
         * These assemblers are called once per cell.
         */
        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> stokes_preconditioner;

        /**
         * A vector of pointers containing all assemblers that compute
         * cell contributions for the Stokes system.
         * These assemblers are called once per cell.
         */
        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> stokes_system;

        /**
         * A vector of pointers containing all assemblers that compute face
         * contributions for the Stokes system. These assemblers are called
         * once per face at a boundary with the properly initialized inputs,
         * therefore they allow terms that only exist on boundary faces (e.g.
         * traction boundary conditions).
         */
        std::vector<std::unique_ptr<Assemblers::Interface<dim>>> stokes_system_on_boundary_face;

        /**
         * A vector of vectors of pointers containing a list of all assemblers
         * for each individual advection system.
         * These assemblers are called once per cell.
         */
        std::vector<std::vector<std::unique_ptr<Assemblers::Interface<dim>>>> advection_system;

        /**
         * A vector of vectors of pointers containing a list of all assemblers
         * for the individual advection systems that compute face contributions
         * at boundaries. These assemblers are called once per boundary face with
         * the properly initialized inputs, therefore they allow terms that only
         * exist on boundary faces (e.g. flux boundary conditions).
         */
        std::vector<std::vector<std::unique_ptr<Assemblers::Interface<dim>>>> advection_system_on_boundary_face;

        /**
         * A vector of vectors of pointers containing a list of all assemblers
         * for the individual advection systems that compute face contributions
         * on faces between cells. These assemblers are called once per interior
         * face with the properly initialized inputs, therefore they allow terms
         * that only exist on interior faces (e.g. DG penalty terms).
         */
        std::vector<std::vector<std::unique_ptr<Assemblers::Interface<dim>>>> advection_system_on_interior_face;

        /**
         * A structure that describes what information an assembler function
         * (listed as one of the assembler objects above) may need to operate.
         *
         * There are a number of pieces of information that are always
         * assumed to be needed. For example, the Stokes and advection
         * assemblers will always need to have access to the material
         * model outputs. But the Stokes assembler may or may not need
         * access to material model outputs for quadrature points on faces.
         *
         * These properties are all preset in a conservative way
         * (i.e., disabled) in the constructor of this class, but can
         * be enabled in Simulator::set_assemblers() when adding
         * individual assemblers. Functions such as
         * Simulator::local_assemble_stokes_preconditioner(),
         * Simulator::local_assemble_stokes_system() will then query
         * these flags to determine whether something has to be
         * initialized for at least one of the assemblers they call.
         */
        struct Properties
        {
          /**
           * Constructor. Disable all properties as described in the
           * class documentation.
           */
          Properties ();

          /**
           * Whether or not at least one of the active assembler objects for
           * a certain equation requires the initialization and re-computation
           * of a MaterialModelOutputs object for each face. This
           * property is only relevant to assemblers that operate on
           * faces.
           */
          bool need_face_material_model_data;

          /**
           * Whether or not at least one of the active assembler objects for
           * a certain equation requires the evaluation of the FEFaceValues
           * object. This is different from need_face_material_model_data,
           * because an assembler might assemble terms that do not require
           * material model outputs.
           */
          bool need_face_finite_element_evaluation;

          /**
           * Whether or not at least one of the active assembler objects for
           * a certain equation requires the computation of the viscosity.
           */
          bool need_viscosity;

          /**
           * A list of FEValues UpdateFlags that are necessary for
           * a given operation. Assembler objects may add to this list
           * as necessary; it will be initialized with a set of
           * "default" flags that will always be set.
           */
          UpdateFlags needed_update_flags;
        };

        /**
         * Lists of properties for the various equations we want to assemble.
         * These property lists are set in Simulator::set_assemblers()
         * where we add individual functions to the vectors of assembler
         * objects above.
         */
        Properties stokes_preconditioner_assembler_properties;
        Properties stokes_system_assembler_properties;
        Properties stokes_system_assembler_on_boundary_face_properties;
        std::vector<Properties> advection_system_assembler_properties;
        std::vector<Properties> advection_system_assembler_on_face_properties;
    };
  }
}


#endif
