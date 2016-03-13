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


#ifndef __aspect__assembly_h
#define __aspect__assembly_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/fe/fe_values.h>


namespace aspect
{
  using namespace dealii;

  template <int dim>
  class Simulator;

  struct AdvectionField;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        struct StokesPreconditioner
        {
          StokesPreconditioner (const FiniteElement<dim> &finite_element,
                                const Quadrature<dim>    &quadrature,
                                const Mapping<dim>       &mapping,
                                const UpdateFlags         update_flags,
                                const unsigned int        n_compositional_fields);
          StokesPreconditioner (const StokesPreconditioner &data);

          virtual ~StokesPreconditioner ();

          FEValues<dim>               finite_element_values;

          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  phi_p;

          std::vector<double>                  temperature_values;
          std::vector<double>                  pressure_values;
          std::vector<SymmetricTensor<2,dim> > strain_rates;
          std::vector<std::vector<double> >    composition_values;

          /**
           * Material model inputs and outputs computed at the current
           * linearization point.
           */
          MaterialModel::MaterialModelInputs<dim> material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> material_model_outputs;
        };



        // We derive the StokesSystem scratch array from the
        // StokesPreconditioner array. We do this because all the objects that
        // are necessary for the assembly of the preconditioner are also
        // needed for the actual system matrix and right hand side, plus some
        // extra data that we need for the time stepping and traction boundaries
        // on the right hand side.
        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &finite_element,
                        const Mapping<dim>       &mapping,
                        const Quadrature<dim>    &quadrature,
                        const Quadrature<dim-1>  &face_quadrature,
                        const UpdateFlags         update_flags,
                        const UpdateFlags         face_update_flags,
                        const unsigned int        n_compositional_fields);

          StokesSystem (const StokesSystem<dim> &data);

          FEFaceValues<dim>               face_finite_element_values;

          std::vector<Tensor<1,dim> >          phi_u;
          std::vector<SymmetricTensor<2,dim> > grads_phi_u;
          std::vector<double>                  div_phi_u;
          std::vector<Tensor<1,dim> >          velocity_values;

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
        };



        template <int dim>
        struct AdvectionSystem
        {
          AdvectionSystem (const FiniteElement<dim> &finite_element,
                           const FiniteElement<dim> &advection_element,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &quadrature,
                           const Quadrature<dim-1>  &face_quadrature,
                           const unsigned int        n_compositional_fields);
          AdvectionSystem (const AdvectionSystem &data);

          FEValues<dim>               finite_element_values;

          std_cxx11::shared_ptr<FEFaceValues<dim> >           face_finite_element_values;
          std_cxx11::shared_ptr<FEFaceValues<dim> >           neighbor_face_finite_element_values;
          std_cxx11::shared_ptr<FESubfaceValues<dim> >        subface_finite_element_values;

          std::vector<types::global_dof_index>   local_dof_indices;

          /**
           * Variables describing the values and gradients of the
           * shape functions at the quadrature points, as they are
           * used in the advection assembly function. note that the sizes
           * of these arrays are equal to the number of shape functions
           * corresponding to the advected field (and not of the overall
           * field!), and that they are also correspondingly indexed.
           */
          std::vector<double>         phi_field;
          std::vector<Tensor<1,dim> > grad_phi_field;
          std::vector<double>         face_phi_field;
          std::vector<Tensor<1,dim> > face_grad_phi_field;
          std::vector<double>         neighbor_face_phi_field;
          std::vector<Tensor<1,dim> > neighbor_face_grad_phi_field;

          std::vector<Tensor<1,dim> > old_velocity_values;
          std::vector<Tensor<1,dim> > old_old_velocity_values;

          std::vector<double>         old_pressure;
          std::vector<double>         old_old_pressure;
          std::vector<Tensor<1,dim> > old_pressure_gradients;
          std::vector<Tensor<1,dim> > old_old_pressure_gradients;

          std::vector<SymmetricTensor<2,dim> > old_strain_rates;
          std::vector<SymmetricTensor<2,dim> > old_old_strain_rates;

          std::vector<double>         old_temperature_values;
          std::vector<double>         old_old_temperature_values;

          std::vector<double>        *old_field_values;
          std::vector<double>        *old_old_field_values;
          std::vector<Tensor<1,dim> > old_field_grads;
          std::vector<Tensor<1,dim> > old_old_field_grads;
          std::vector<double>         old_field_laplacians;
          std::vector<double>         old_old_field_laplacians;

          std::vector<std::vector<double> > old_composition_values;
          std::vector<std::vector<double> > old_old_composition_values;

          std::vector<double>         current_temperature_values;
          std::vector<Tensor<1,dim> > current_velocity_values;
          std::vector<Tensor<1,dim> > face_current_velocity_values;
          std::vector<Tensor<1,dim> > mesh_velocity_values;
          std::vector<Tensor<1,dim> > face_mesh_velocity_values;

          std::vector<SymmetricTensor<2,dim> > current_strain_rates;
          std::vector<std::vector<double> > current_composition_values;

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
           * Material model inputs and outputs computed at a previous
           * time step's solution, or an extrapolation from previous
           * time steps.
           */
          MaterialModel::MaterialModelInputs<dim> explicit_material_model_inputs;
          MaterialModel::MaterialModelOutputs<dim> explicit_material_model_outputs;

          /**
           * Heating model outputs computed at the quadrature points of the
           * current cell at the time of the current linearization point.
           * As explained in the class documentation of
           * HeatingModel::HeatingModelOutputs each term contains the sum of all
           * enabled heating mechanism contributions.
           */
          HeatingModel::HeatingModelOutputs heating_model_outputs;
        };



      }


      // The CopyData arrays are similar to the
      // Scratch arrays. They provide a
      // constructor, a copy operation, and
      // some arrays for local matrix, local
      // vectors and the relation between local
      // and global degrees of freedom (a.k.a.
      // <code>local_dof_indices</code>).
      namespace CopyData
      {
        template <int dim>
        struct StokesPreconditioner
        {
          StokesPreconditioner (const FiniteElement<dim> &finite_element);
          StokesPreconditioner (const StokesPreconditioner &data);

          virtual ~StokesPreconditioner ();

          FullMatrix<double>          local_matrix;
          std::vector<types::global_dof_index>   local_dof_indices;
        };



        template <int dim>
        struct StokesSystem : public StokesPreconditioner<dim>
        {
          StokesSystem (const FiniteElement<dim> &finite_element,
                        const bool                do_pressure_rhs_compatibility_modification);
          StokesSystem (const StokesSystem<dim> &data);

          Vector<double> local_rhs;
          Vector<double> local_pressure_shape_function_integrals;
        };



        template <int dim>
        struct AdvectionSystem
        {
          /**
           * Constructor.
           * @param finite_element The element that describes the field for which we
           *    are trying to assemble a linear system. <b>Not</b> the global finite
           *    element.
           */
          AdvectionSystem (const FiniteElement<dim> &finite_element,
                           const bool                field_is_discontinuous);
          AdvectionSystem (const AdvectionSystem &data);

          /**
           * Local contributions to the global matrix and right hand side
           * that correspond only to the variables listed in local_dof_indices
           */
          FullMatrix<double>          local_matrix;
          /** Local contributions to the global matrix from the face terms in the
           * discontinuous Galerkin method. The vectors are of length
           * GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell
           * so as to hold one matrix for each possible face or subface of the cell.
           * The discontinuous Galerkin bilinear form contains terms arising from internal
           * (to the cell) values and external (to the cell) values.
           * _int_ext and ext_int hold the terms arising from the pairing between a cell
           * and its neighbor, while _ext_ext is the pairing of the neighbor's dofs with
           * themselves. In the continuous Galerkin case, these are unused, and set to size zero.
           **/
          std::vector<FullMatrix<double> >         local_matrices_int_ext;
          std::vector<FullMatrix<double> >         local_matrices_ext_int;
          std::vector<FullMatrix<double> >         local_matrices_ext_ext;
          Vector<double>              local_rhs;

          /** Denotes which face matrices have actually been assembled in the DG field
           * assembly. Entries for matrices not used (for example, those corresponding
           * to non-existent subfaces; or faces being assembled by the neighboring cell)
           * are set to false.
           **/
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
          /** Indices of the degrees of freedom corresponding to the temperature
           * or composition field on all possible neighboring cells. This is used
           * in the discontinuous Galerkin method. The outer std::vector has
           *  length GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell,
           * and has size zero if in the continuous Galerkin case.
           **/
          std::vector<std::vector<types::global_dof_index> >   neighbor_dof_indices;
        };
      }

      /**
       * A namespace for the definition of sub-namespaces in which we
       * implement assemblers for the various terms in the linear systems
       * ASPECT solves.
       */
      namespace Assemblers
      {
        /**
         * A base class for objects that implement assembly
         * operations.
         *
         * The point of this class is primarily so that we can store
         * pointers to such objects in a list. The objects are created
         * in Simulator::set_assemblers() (which is the only place where
         * we know their actual type) and destroyed in the destructor of
         * class Simulation.
         */
        template <int dim>
        class AssemblerBase
        {
          public:
            virtual ~AssemblerBase () {}

            /**
             * This function gets called if a MaterialModelOutputs is created
             * and allows the assembler to attach AdditionalOutputs. The
             * function might be called more than once for a
             * MaterialModelOutput, so it is recommended to check if
             * get_additional_output() returns an instance before adding a new
             * one to the additional_outputs vector.
             */
            virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &) {}
        };
      }

      /**
       * A structure that consists of member variables representing
       * each one or more functions that need to be called when
       * assembling right hand side vectors, matrices, or complete linear
       * systems. We use this approach in order to support the following
       * cases:
       * - Assembling different formulations: When assembling either the
       *   full equations or only the Boussinesq approximation (to give just
       *   two examples), one needs different terms. This could be achieved
       *   using a large number of <code>switch</code> or <code>if</code>
       *   statements in the code, or one could encapsulate each equation
       *   or approximation into a collection of functions for this particular
       *   purpose. The approach chosen here in essence allows the
       *   implementation of each set of equations in its own scope, and we
       *   then just need to store a pointer to the function that
       *   assembles the Stokes system (for example) for the selected
       *   approximation. The pointer to this function is stored in the
       *   appropriate member variable of this class.
       * - Sometimes, we want to assemble a number of terms that build on
       *   each other. An example is the addition of free boundary terms
       *   to the Stokes matrix. Rather than having to "know" in one
       *   place about all of the terms that need to be assembled,
       *   we simply add the function that computes these terms as
       *   another "slot" to the appropriate "signal" declared in this
       *   class.
       */
      template <int dim>
      struct AssemblerLists
      {
        /**
         * A structure used to accumulate the results of the
         * compute_advection_system_residual slot functions below.
         * It takes an iterator range of vectors and returns the element-wise
         * sum of the values.
         */
        template<typename VectorType>
        struct ResidualWeightSum
        {
          typedef VectorType result_type;

          template<typename InputIterator>
          VectorType operator()(InputIterator first, InputIterator last) const
          {
            // If there are no slots to call, just return the
            // default-constructed value
            if (first == last)
              return VectorType();

            VectorType sum = *first++;
            while (first != last)
              {
                Assert(sum.size() == first->size(),
                       ExcMessage("Invalid vector length for residual summation."));

                for (unsigned int i = 0; i < first->size(); ++i)
                  sum[i] += (*first)[i];

                ++first;
              }

            return sum;
          }
        };


        /**
         * A signal that is called from Simulator::local_assemble_stokes_preconditioner()
         * and whose slots are supposed to assemble terms that together form the
         * Stokes preconditioner matrix.
         *
         * The arguments to the slots are as follows:
         * - The Simulator::pressure_scaling value used to scale velocity
         *   and pressure components against each other.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const double,
                                      internal::Assembly::Scratch::StokesPreconditioner<dim>  &,
                                      internal::Assembly::CopyData::StokesPreconditioner<dim> &)> local_assemble_stokes_preconditioner;

        /**
         * A signal that is called from Simulator::local_assemble_stokes_system()
         * and whose slots are supposed to assemble terms that together form the
         * Stokes system matrix and right hand side.
         *
         * The arguments to the slots are as follows:
         * - The cell on which we currently assemble.
         * - The Simulator::pressure_scaling value used to scale velocity
         *   and pressure components against each other.
         * - Whether or not to actually assemble the matrix. If @p false,
         *   then only assemble the right hand side of the Stokes system.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const typename DoFHandler<dim>::active_cell_iterator &,
                                      const double,
                                      const bool,
                                      internal::Assembly::Scratch::StokesSystem<dim>       &,
                                      internal::Assembly::CopyData::StokesSystem<dim>      &)> local_assemble_stokes_system;

        /**
         * A signal that is called from Simulator::local_assemble_stokes_system()
         * and whose slots are supposed to assemble terms that together form the
         * Stokes system matrix and right hand side. This signal is called
         * once for each boundary face.
         *
         * The arguments to the slots are as follows:
         * - The cell on which we currently assemble.
         * - The number of the face on which we intend to assemble. This
         *   face (of the current cell) will be at the boundary of the
         *   domain.
         * - The Simulator::pressure_scaling value used to scale velocity
         *   and pressure components against each other.
         * - Whether or not to actually assemble the matrix. If @p false,
         *   then only assemble the right hand side of the Stokes system.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const typename DoFHandler<dim>::active_cell_iterator &,
                                      const unsigned int,
                                      const double,
                                      const bool,
                                      internal::Assembly::Scratch::StokesSystem<dim>       &,
                                      internal::Assembly::CopyData::StokesSystem<dim>      &)> local_assemble_stokes_system_on_boundary_face;

        /**
         * A signal that is called from Simulator::local_assemble_advection_system()
         * and whose slots are supposed to assemble terms that together form the
         * Advection system matrix and right hand side.
         *
         * The arguments to the slots are as follows:
         * - The cell on which we currently assemble.
         * - The advection field that is currently assembled.
         * - The artificial viscosity computed for the current cell
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * - The copy object into which assemblers add up their contributions.
         */
        boost::signals2::signal<void (const typename DoFHandler<dim>::active_cell_iterator &,
                                      const typename Simulator<dim>::AdvectionField &,
                                      const double,
                                      internal::Assembly::Scratch::AdvectionSystem<dim>       &,
                                      internal::Assembly::CopyData::AdvectionSystem<dim>      &)> local_assemble_advection_system;

        /**
         * A signal that is called from Simulator::local_assemble_advection_system()
         * and whose slots are supposed to compute an advection system residual
         * according to the terms their assemblers implement. The residuals are
         * computed and returned in a std::vector<double> with one entry per
         * quadrature point of this cell. If more than one
         * assembler is connected the residuals are simply added up by the
         * ResidualWeightSum structure.
         *
         * The arguments to the slots are as follows:
         * - The cell on which we currently assemble.
         * - The advection field that is currently assembled.
         * - The scratch object in which temporary data is stored that
         *   assemblers may need.
         * The return value of each slot is:
         * - A std::vector<double> containing one residual entry per quadrature
         *   point.
         */
        boost::signals2::signal<std::vector<double> (const typename DoFHandler<dim>::active_cell_iterator &,
                                                     const typename Simulator<dim>::AdvectionField &,
                                                     internal::Assembly::Scratch::AdvectionSystem<dim> &),
                                                              ResidualWeightSum<std::vector<double> > > compute_advection_system_residual;

        /**
         * A structure that describes what information an assembler function
         * (listed as one of the signals/slots above) may need to operate.
         *
         * There are a number of pieces of information that are always
         * assumed to be needed. For example, the Stokes and advection
         * assemblers will always need to have access to the material
         * model outputs. But the Stokes assembler may or may not need
         * access to material model output for quadrature points on faces.
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
           * Whether or not at least one of the the assembler slots in
           * a signal require the initialization and re-computation of
           * a MaterialModelOutputs object for each face. This
           * property is only relevant to assemblers that operate on
           * boundary faces.
           */
          bool need_face_material_model_data;

          /**
           * A list of FEValues UpdateFlags that are necessary for
           * a given operation. Assembler objects may add to this list
           * as necessary; it will be initialized with a set of
           * "default" flags that will always be set.
           */
          UpdateFlags needed_update_flags;
        };

        /**
         * A list of properties of the various types of assemblers.
         * These property lists are set in Simulator::set_assemblers()
         * where we add individual functions to the signals above.
         */
        Properties stokes_preconditioner_assembler_properties;
        Properties stokes_system_assembler_properties;
        Properties stokes_system_assembler_on_boundary_face_properties;
      };

    }
  }
}

#endif
