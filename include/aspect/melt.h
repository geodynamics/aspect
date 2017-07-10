/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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
#include <aspect/assembly.h>
#include <aspect/material_model/interface.h>
#include <aspect/boundary_fluid_pressure/interface.h>

namespace aspect
{
  using namespace dealii;

  namespace MaterialModel
  {
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
                      const FullMatrix<double>  &expansion_matrix);
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
    };


  }



  namespace Assemblers
  {
    /**
     * A class for the definition of functions that implement the
     * linear system terms for the *melt* migration compressible or
     * incompressible equations.
     */
    template <int dim>
    class MeltEquations : public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Attach melt outputs.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;


        /**
         * Compute the integrals for the preconditioner for the Stokes system in
         * the case of melt migration on a single cell.
         */
        void
        local_assemble_stokes_preconditioner_melt (const double                                             pressure_scaling,
                                                   internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                                                   internal::Assembly::CopyData::StokesPreconditioner<dim> &data) const;

        /**
         * Compute the integrals for the Stokes matrix and right hand side in
         * the case of melt migration on a single cell.
         */
        void
        local_assemble_stokes_system_melt (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           const double                                     pressure_scaling,
                                           const bool                                       rebuild_stokes_matrix,
                                           internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                           internal::Assembly::CopyData::StokesSystem<dim> &data) const;

        /**
         * Compute the boundary integrals for the Stokes right hand side in
         * the case of melt migration on a single cell. These boundary terms
         * are used to describe Neumann boundary conditions for the fluid pressure.
         */
        void
        local_assemble_stokes_system_melt_boundary (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                    const unsigned int                                    face_no,
                                                    const double                                          pressure_scaling,
                                                    internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                                    internal::Assembly::CopyData::StokesSystem<dim>      &data) const;

        /**
         * Compute the integrals for the Advection system matrix and right hand side in
         * the case of melt migration on a single cell.
         */
        void
        local_assemble_advection_system_melt (const typename Simulator<dim>::AdvectionField      &advection_field,
                                              const double                                        artificial_viscosity,
                                              internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                                              internal::Assembly::CopyData::AdvectionSystem<dim> &data) const;

        /**
         * Compute the residual of the advection system on a single cell in
         * the case of melt migration.
         */
        std::vector<double>
        compute_advection_system_residual_melt(const typename Simulator<dim>::AdvectionField     &advection_field,
                                               internal::Assembly::Scratch::AdvectionSystem<dim> &scratch) const;


      private:
        /**
         * Compute the right-hand side of the fluid pressure equation of the Stokes
         * system in case the simulation uses melt transport. This includes a term
         * derived from Darcy's law, a term including the melting rate and a term dependent
         * on the densities and velocities of fluid and solid.
         */
        double
        compute_fluid_pressure_rhs(const internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                   MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                   MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                                   const unsigned int q_point) const;

        /**
         * Compute the right-hand side for the porosity system. This function
         * returns 0 for the temperature and all of the compositional fields, except
         * for the porosity. It includes the melting rate and a term dependent
         * on the density and velocity.
         *
         * This function is implemented in
         * <code>source/simulator/melt.cc</code>.
         */
        double
        compute_melting_RHS(const internal::Assembly::Scratch::AdvectionSystem<dim>  &scratch,
                            const typename MaterialModel::Interface<dim>::MaterialModelInputs &material_model_inputs,
                            const typename MaterialModel::Interface<dim>::MaterialModelOutputs &material_model_outputs,
                            const typename Simulator<dim>::AdvectionField &advection_field,
                            const unsigned int q_point) const;
    };

    /**
     * A namespace for the definition of functions that implement various
     * other terms that need to occasionally or always be assembled.
     */
    namespace OtherTerms
    {
      /**
       * Integrate the local fluid pressure shape functions on a single cell
       * for models with melt migration, so that they can later be used to do
       * the pressure right-hand side compatibility modification.
       */
      template <int dim>
      void
      pressure_rhs_compatibility_modification_melt (const SimulatorAccess<dim>                      &simulator_access,
                                                    internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                    internal::Assembly::CopyData::StokesSystem<dim> &data);
    }
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
       * Declare additional parameters that are needed in models with
       * melt transport (including the fluid pressure boundary conditions).
       */
      static void declare_parameters (ParameterHandler &prm);

      /**
       * Parse additional parameters that are needed in models with
       * melt transport (including the fluid pressure boundary conditions).
       *
       * This has to be called before edit_finite_element_variables,
       * so that the finite elements that are used for the additional melt
       * variables can be specified in the input file and are parsed before
       * the introspection object is created.
       */
      void parse_parameters (ParameterHandler &prm);

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
       * Setup SimulatorAccess for the plugins related to melt transport.
       */
      void initialize_simulator (const Simulator<dim> &simulator_object);

      /**
       * Compute fluid velocity and solid pressure in this ghosted solution vector.
       * The fluid velocity is computed by solving a mass matrix problem, and the
       * solid pressure is computed algebraically.
       *
       * @param solution The existing solution vector that contains the values
       * for porosity, compaction pressure, fluid pressure and solid velocity
       * obtained by solving the Stokes and advection system, and that will be
       * updated with the computed values for fluid velocity and solid pressure.
       */
      void compute_melt_variables(LinearAlgebra::BlockVector &solution);

      /**
       * Return whether this object refers to the porosity field.
       */
      bool is_porosity (const typename Simulator<dim>::AdvectionField &advection_field) const;

      /**
       * The porosity limit for melt migration. For smaller porosities, the equations
       * reduce to the Stokes equations and neglect melt transport. In practice, this
       * means that all terms in the assembly related to the migration of melt are set
       * to zero for porosities smaller than this threshold.
       * This does not include the compaction term $p_c/\xi$, which is necessary for the
       * solvability of the linear system, but does not influence the solution variables
       * of the Stokes problem (in the absence of porosity).
       */
      double melt_transport_threshold;

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
       * Store a pointer to the fluid pressure boundary plugin, so that the
       * initialization can be done together with the other objects related to melt
       * transport.
       */
      std_cxx11::unique_ptr<BoundaryFluidPressure::Interface<dim> > boundary_fluid_pressure;
  };

}

#endif
