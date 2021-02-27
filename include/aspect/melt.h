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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
 */


#ifndef _aspect_melt_h
#define _aspect_melt_h

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/simulator/assemblers/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>

namespace aspect
{
  using namespace dealii;

  namespace MaterialModel
  {
    /**
     * The MeltInputs provide the compaction pressures and
     * melt (fluid) velocities, so that they can be used as
     * additional inputs in heating or material models.
     */
    template <int dim>
    class MeltInputs : public AdditionalMaterialInputs<dim>
    {
      public:
        /**
         * Constructor. When the MeltInputs are created,
         * all properties are initialized with signalingNaNs.
         * This means that individual heating or material models
         * can all attach the plugins they need, and in a later
         * step they will all be filled together (using the fill
         * function).
         */
        MeltInputs (const unsigned int n_points);

        /**
         * Compaction pressure values $p_c$ at the given positions.
         */
        std::vector<double> compaction_pressures;

        /**
         * An approximation for the fluid (melt) velocities
         * at the given positions.
         */
        std::vector<Tensor<1,dim> > fluid_velocities;

        /**
         * Fill the compaction pressures and fluid velocities.
         */
        void fill (const LinearAlgebra::BlockVector &solution,
                   const FEValuesBase<dim>          &fe_values,
                   const Introspection<dim>         &introspection) override;
    };

    template <int dim>
    class MeltOutputs : public AdditionalMaterialOutputs<dim>
    {
      public:
        MeltOutputs (const unsigned int n_points,
                     const unsigned int /*n_comp*/)
        {
          compaction_viscosities.resize(n_points);
          fluid_viscosities.resize(n_points);
          permeabilities.resize(n_points);
          fluid_densities.resize(n_points);
          fluid_density_gradients.resize(n_points, Tensor<1,dim>());
        }

        /**
         * Compaction viscosity values $\xi$ at the given positions.
         * This parameter describes the resistance of the solid matrix
         * in a two-phase simulation to dilation and compaction.
         */
        std::vector<double> compaction_viscosities;

        /**
         * Fluid (melt) viscosity values $\eta_f$ at the given positions.
         */
        std::vector<double> fluid_viscosities;

        /**
         * Permeability values $k$ at the given positions.
         */
        std::vector<double> permeabilities;

        /**
         * Fluid (melt) density values $\rho_f$ at the given positions.
         */
        std::vector<double> fluid_densities;

        /**
         * An approximation for the fluid (melt) density gradients
         * $\nabla \rho_f$ at the given positions. These values are
         * required for compressible models to describe volume changes
         * of melt in dependence of pressure, temperature etc.
         */
        std::vector<Tensor<1,dim> > fluid_density_gradients;

        /**
         * Do the requested averaging operation for the melt outputs.
         * The projection matrix argument is only used if the operation
         * chosen is project_to_Q1.
         */
        void average (const MaterialAveraging::AveragingOperation operation,
                      const FullMatrix<double>  &projection_matrix,
                      const FullMatrix<double>  &expansion_matrix) override;
    };

    /**
     * Base class for material models that implement a melt fraction function.
     * This is used to compute some statistics about the melt fraction.
     */
    template <int dim>
    class MeltFractionModel
    {
      public:
        /**
         * Compute the equilibrium melt fractions for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * equilibrium melt fraction for each given input conditions.
         */
        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const = 0;

        /**
         * Destructor. Does nothing but is virtual so that derived classes
         * destructors are also virtual.
         */
        virtual ~MeltFractionModel ()
        {};
    };

    /**
     * Base class for material models to be used with melt transport enabled.
     */
    template <int dim>
    class MeltInterface: public MaterialModel::Interface<dim>
    {
      public:
        /**
         * Reference value for the Darcy coefficient, which is defined as
         * permeability divided by fluid viscosity. Units: m^2/Pa/s.
         */
        virtual double reference_darcy_coefficient () const = 0;

        /**
         * Returns the cell-averaged and cut-off value of p_c_scale,
         * the factor we use to rescale the compaction pressure and to
         * decide if a cell is a melt cell.
         * The last input argument @p consider_is_melt_cell determines if
         * this computation takes into account if a cell is a "melt cell".
         * Melt cells are cells where we solve the melt transport equations,
         * as indicated by the entries stored in the is_melt_cell vector of
         * the melt handler. In case @p consider_is_melt_cell is set to true,
         * this function returns a value of zero if the cell is not a melt cell.
         * If @p consider_is_melt_cell is set to false the computation
         * disregards the information about which cells are melt cells,
         * and computes p_c_scale from the cell-averaged Darcy coefficient
         * for all cells. This is needed for example when we want to update
         * the is_melt_cell vector and need to find out which cells should be
         * marked as melt cells.
         */
        double p_c_scale (const MaterialModel::MaterialModelInputs<dim> &inputs,
                          const MaterialModel::MaterialModelOutputs<dim> &outputs,
                          const MeltHandler<dim> &melt_handler,
                          const bool consider_is_melt_cell) const;
    };


  }



  namespace Assemblers
  {
    /**
     * A base class for the definition of assemblers that implement the linear
     * system terms for the *melt* migration compressible or incompressible
     * equations.
     */
    template <int dim>
    class MeltInterface : public aspect::Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Attach melt outputs. Since most melt assemblers require the
         * melt material model properties they are created in this base class
         * already.
         */
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * Compute the integrals for the preconditioner for the Stokes system in
     * the case of melt migration on a single cell.
     */
    template <int dim>
    class MeltStokesPreconditioner : public MeltInterface<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * Compute the integrals for the Stokes matrix and right hand side in
     * the case of melt migration on a single cell.
     */
    template <int dim>
    class MeltStokesSystem : public MeltInterface<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };


    /**
     * Compute the boundary integrals for the Stokes right hand side in
     * the case of melt migration on a single cell. These boundary terms
     * are used to describe Neumann boundary conditions for the fluid pressure.
     */
    template <int dim>
    class MeltStokesSystemBoundary : public MeltInterface<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * Compute the integrals for the Advection system matrix and right hand side in
     * the case of melt migration on a single cell.
     */
    template <int dim>
    class MeltAdvectionSystem : public MeltInterface<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        /**
         * Compute the residual of the advection system on a single cell in
         * the case of melt migration.
         */
        std::vector<double>
        compute_residual(internal::Assembly::Scratch::ScratchBase<dim> &scratch_base) const override;
    };

    /**
     * Integrate the local fluid pressure shape functions on a single cell
     * for models with melt migration, so that they can later be used to do
     * the pressure right-hand side compatibility modification.
     */
    template <int dim>
    class MeltPressureRHSCompatibilityModification : public MeltInterface<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * Assemble traction boundary condition terms for models with melt.
     */
    template <int dim>
    class MeltBoundaryTraction : public MeltInterface<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };
  }


  namespace Melt
  {
    template <int dim>
    struct Parameters
    {
      /**
       * Declare additional parameters that are needed in models with
       * melt transport.
       */
      static void declare_parameters (ParameterHandler &prm);

      /**
       * Parse additional parameters that are needed in models with
       * melt transport.
       *
       * This has to be called before edit_finite_element_variables,
       * so that the finite elements that are used for the additional melt
       * variables can be specified in the input file and are parsed before
       * the introspection object is created.
       */
      void parse_parameters (ParameterHandler &prm);

      /**
       * The factor by how much the Darcy coefficient K_D in a cell can be smaller than
       * the reference Darcy coefficient for this cell still to be considered a melt cell
       * (for which the melt transport equations are solved). If the Darcy coefficient
       * is smaller than the product of this value and the reference Dracy coefficient,
       * the cell is not considered a melt cell and the Stokes system without melt
       * transport is solved instead. In practice, this means that all terms in the
       * assembly related to the migration of melt are set to zero, and the compaction
       * pressure degrees of freedom are constrained to be zero.
       */
      double melt_scaling_factor_threshold;

      /**
       * Whether to use a porosity weighted average of the melt and solid velocity
       * to advect heat in the temperature equation or not. If this is set to true,
       * additional terms are assembled on the left-hand side of the temperature
       * advection equation in models with melt migration.
       * If this is set to false, only the solid velocity is used (as in models
       * without melt migration).
       */
      bool heat_advection_by_melt;

      /**
       * Whether to use a discontinuous element for the compaction pressure or not.
       */
      bool use_discontinuous_p_c;

      /**
       * Whether to cell-wise average the material properties that are used to
       * compute the melt velocity or not. Note that the melt velocity is computed
       * as the sum of the solid velocity and the phase separation flux (difference
       * between melt and solid velocity). If this parameter is set to true,
       * material properties in the computation of the phase separation flux will
       * be averaged cell-wise.
       */
      bool average_melt_velocity;
    };
  }


  /**
   * Class that contains all runtime parameters and other helper functions
   * related to melt transport. A global instance can be retrieved with
   * SimulatorAccess<dim>::get_melt_handler(), but keep in mind that it only
   * exists if parameters.include_melt_transport is true.
   */
  template <int dim>
  class MeltHandler: public SimulatorAccess<dim>
  {
    public:
      MeltHandler(ParameterHandler &prm);

      /**
       * Create an additional material model output object that contains
       * the additional output variables needed in simulation with melt transport,
       * and attaches a pointer to it to the corresponding vector in the
       * MaterialModel::MaterialModelOutputs structure.
       */
      static void create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output);

      /**
       * Add the additional variables we need in simulations with melt
       * migration to the list of variables, which will be used later
       * to set up the introspection object.
       */
      void edit_finite_element_variables(const Parameters<dim> &parameters,
                                         std::vector<VariableDeclaration<dim> > &variables);

      /**
       * Determine, based on the run-time parameters of the current simulation,
       * which functions need to be called in order to assemble linear systems,
       * matrices, and right hand side vectors.
       */
      void set_assemblers (Assemblers::Manager<dim> &assemblers) const;


      /**
       * Initialize function. This is mainly to check that the melt transport
       * parameters chosen in the input file are consistent with the rest of
       * the options. We can not do this in the parse_parameters function,
       * as we do not have simulator access at that point.
       */
      void initialize() const;

      /**
       * Setup SimulatorAccess for the plugins related to melt transport.
       */
      void initialize_simulator (const Simulator<dim> &simulator_object) override;

      /**
       * Compute fluid velocity and solid pressure in this ghosted solution vector.
       * The fluid velocity is computed by solving a mass matrix problem, and the
       * solid pressure is computed algebraically.
       *
       * @param system_matrix The system matrix with an already set up sparsity
       * pattern that will be used by this function to compute the melt variables.
       * @param solution The existing solution vector that contains the values
       * for porosity, compaction pressure, fluid pressure and solid velocity
       * obtained by solving the Stokes and advection system, and that will be
       * updated with the computed values for fluid velocity and solid pressure.
       * @param system_rhs The right-hand side vector that will be used by
       * this function to compute the melt variables.
       */
      void compute_melt_variables(LinearAlgebra::BlockSparseMatrix &system_matrix,
                                  LinearAlgebra::BlockVector &solution,
                                  LinearAlgebra::BlockVector &system_rhs);

      /**
       * Return whether this object refers to the porosity field.
       */
      bool is_porosity (const typename Simulator<dim>::AdvectionField &advection_field) const;

      /**
       * Apply free surface stabilization to a cell of the system matrix when melt
       * transport is used in the computation. Called during assembly of the system matrix.
       */
      void apply_free_surface_stabilization_with_melt (const double free_surface_theta,
                                                       const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                       internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                                       internal::Assembly::CopyData::StokesSystem<dim>      &data) const;

      /**
       * Constrain the compaction pressure to zero in all cells that are not
       * "melt cells" (cells where the porosity is above a given threshold).
       * This reverts the system of equations we solve back to the Stokes
       * system without melt transport for these cells.
       */
      void add_current_constraints(AffineConstraints<double> &constraints);

      /**
       * Returns the entry of the private variable is_melt_cell_vector for the
       * cell given in the input, describing if we have melt transport in this
       * cell or not.
       */
      bool is_melt_cell(const typename DoFHandler<dim>::active_cell_iterator &cell) const;

      /**
       * Given the Darcy coefficient @p K_D as computed by the material model,
       * limit the coefficient to a minimum value (computed as the K_D
       * variation threshold given in the input file times the reference Darcy
       * coefficient) in melt cells and return this value. If @p is_melt_cell
       * is false, return zero.
       */
      double limited_darcy_coefficient(const double K_D,
                                       const bool is_melt_cell) const;

      /**
       * Return a pointer to the boundary fluid pressure.
       */
      const BoundaryFluidPressure::Interface<dim> &
      get_boundary_fluid_pressure () const;

      /**
       * The object that stores the run-time parameters that control the how the
       * melt transport equations are solved.
       */
      Melt::Parameters<dim> melt_parameters;

    private:
      /**
       * Store a pointer to the fluid pressure boundary plugin, so that the
       * initialization can be done together with the other objects related to melt
       * transport.
       */
      const std::unique_ptr<aspect::BoundaryFluidPressure::Interface<dim> > boundary_fluid_pressure;

      /**
       * is_melt_cell_vector[cell->active_cell_index()] says whether we want to
       * solve the melt transport equations (as opposed to the Stokes equations without
       * melt) in this cell or not. The value is set to true or false based on the
       * porosity in the cell (only cells where the porosity is above a threshold are
       * considered melt cells).
       */
      std::vector<bool> is_melt_cell_vector;

      /**
       * Constraint object. We need to save the current constraints at the
       * start of every time step so that we can add the melt constraints,
       * which depend on the solution of the porosity field, later after
       * we have computed this solution.
       */
      AffineConstraints<double> current_constraints;

  };

}

#endif
