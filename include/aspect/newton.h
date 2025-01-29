/*
  Copyright (C) 2016 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_newton_h
#define _aspect_newton_h

#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/material_model/interface.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * This class holds the derivatives for the Newton solver.
     */
    template <int dim>
    class MaterialModelDerivatives : public AdditionalMaterialOutputs<dim>
    {
      public:
        /**
         * Constructor. Initialize the various arrays of this structure with the
         * given number of quadrature points.
         */
        MaterialModelDerivatives (const unsigned int n_points);

        /**
         * The derivatives of the viscosities with respect to pressure.
         */
        std::vector<double> viscosity_derivative_wrt_pressure;

        /**
         * The derivatives of the viscosities with respect to strain rate.
         */
        std::vector<SymmetricTensor<2,dim>> viscosity_derivative_wrt_strain_rate;

        /**
         * The weights used for calculating the averages of viscosity
         * derivatives when material averaging is applied.
         */
        std::vector<double> viscosity_derivative_averaging_weights;
    };
  }


  namespace Newton
  {
    struct Parameters
    {

      /**
       * This enum describes the type of stabilization is used
       * for the Newton solver. None represents no stabilization,
       * SPD represent that the resulting matrix is made Symmetric
       * Positive Definite, symmetric represents that the matrix is
       * only symmetrized, and PD represents that we do the same as
       * what we do for SPD, but without the symmetrization.
       */
      enum Stabilization
      {
        none = 0,
        symmetric = 1,
        PD = 2,
        SPD = symmetric | PD
      };


      /**
       * A binary 'or' operator that concatenates a set of stabilization
       * flags by returning an object that combines the bits set in each
       * of the two arguments.
       */
      friend
      Stabilization
      operator| (const Stabilization a,
                 const Stabilization b)
      {
        return static_cast<Stabilization>(
                 static_cast<int>(a) | static_cast<int>(b));
      }


      /**
       * A binary 'and' operator that takes the intersection of two sets
       * of stabilization flags by returning an object that selects those bits
       * that are set in both of the two arguments.
       */
      friend
      Stabilization
      operator& (const Stabilization a,
                 const Stabilization b)
      {
        return static_cast<Stabilization>(
                 static_cast<int>(a) & static_cast<int>(b));
      }


      /**
       * Declare additional parameters that are needed for the Newton.
       * solver.
       */
      static void declare_parameters (ParameterHandler &prm);

      /**
       * Parse additional parameters that are needed for the Newton.
       * solver.
       */
      void parse_parameters (ParameterHandler &prm);

      /**
       * A scaling factor used for scaling the
       * derivative part of the Newton Stokes solver in the assembly.
       *
       * The exact Newton matrix consists of the Stokes matrix plus a term
       * that results from the linearization of the material coefficients.
       * The scaling factor multiplies these additional terms. In a full
       * Newton method, it would be equal to one, but it can be chosen
       * smaller in cases where the resulting linear system has undesirable
       * properties.
       *
       * If the scaling factor is zero, the resulting matrix is simply the
       * Stokes matrix, and the resulting scheme is a defect correction
       * (i.e., Picard iteration).
       */
      double              newton_derivative_scaling_factor;

      Stabilization       preconditioner_stabilization;
      Stabilization       velocity_block_stabilization;

      /**
       * Whether to use the Newton failsafe or not. If the failsafe is used, a failure
       * of the linear solver is caught and we try to solve it again with both the
       * preconditioner and the velocity block being stabilized with the SPD stabilization.
       */
      bool                use_Newton_failsafe;

      /**
       * The nonlinear tolerance at which to switch the
       * nonlinear solver from defect correction Picard to
       * Newton.
       */
      double              nonlinear_switch_tolerance;

      bool                use_Eisenstat_Walker_method_for_Picard_iterations;
      unsigned int        max_pre_newton_nonlinear_iterations;
      unsigned int        max_newton_line_search_iterations;
      bool                use_newton_residual_scaling_method;
      double              maximum_linear_stokes_solver_tolerance;
      double              SPD_safety_factor;
    };


    /**
     * Get a std::string describing the stabilization type used for the
     * preconditioner.
     */
    std::string
    to_string(const Newton::Parameters::Stabilization preconditioner_stabilization);
  }



  /**
   * A class that supports the functionality of the Newton solver.
   */
  template <int dim>
  class NewtonHandler : public SimulatorAccess<dim>
  {
    public:
      /**
       * Determine, based on the run-time parameters of the current simulation,
       * which functions need to be called in order to assemble linear systems,
       * matrices, and right hand side vectors.
       */
      void set_assemblers (Assemblers::Manager<dim> &assemblers) const;

      /**
       * Create an additional material model output object that contains
       * the additional output variables (the derivatives) needed for the
       * Newton solver.
       */
      static void create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output);

      /**
       * The object that stores the run-time parameters that control the Newton
       * method.
       */
      Newton::Parameters parameters;
  };


  namespace Assemblers
  {
    /**
     * A base class for the definition of assemblers that implement the linear
     * system terms for the NewtonStokes solver scheme.
     */
    template <int dim>
    class NewtonInterface : public aspect::Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Attach Newton outputs. Since most Newton assemblers require the
         * material model derivatives they are created in this base class
         * already.
         */
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * This class assembles the terms of the Newton Stokes preconditioner matrix for the current cell.
     */
    template <int dim>
    class NewtonStokesPreconditioner : public NewtonInterface<dim>
    {
      public:
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        /**
         * Create additional material models outputs for computing viscoelastic strain rate when
         * elasticity is enabled.
         */
        void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Newton Stokes system for the current cell.
     */
    template <int dim>
    class NewtonStokesIncompressibleTerms : public NewtonInterface<dim>
    {
      public:
        void
        execute (internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        /**
         * Create additional material models outputs for assembly of derivatives or adding additional
         * terms to the right hand side of the Stokes equations. The latter could include viscoelastic
         * forces or other user-defined values calculated within the material model.
         */
        void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * This class assembles the term that arises in the viscosity term of the Newton Stokes matrix for
     * compressible models, because the divergence of the velocity is not longer zero.
     */
    template <int dim>
    class NewtonStokesCompressibleStrainRateViscosityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the right-hand-side term of the Newton Stokes system
     * that is caused by the compressibility in the mass conservation equation.
     * This function approximates this term as
     * $- \nabla \mathbf{u} = \frac{1}{\rho} \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u}$
     */
    template <int dim>
    class NewtonStokesReferenceDensityCompressibilityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the compressibility term of the Newton Stokes system
     * that is caused by the compressibility in the mass conservation equation.
     * It includes this term implicitly in the matrix,
     * which is therefore not longer symmetric.
     * This function approximates this term as
     * $ - \nabla \mathbf{u} - \frac{1}{\rho} \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u} = 0$
     */
    template <int dim>
    class NewtonStokesImplicitReferenceDensityCompressibilityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the right-hand-side term of the Newton Stokes system
     * that is caused by the compressibility in the mass conservation equation.
     * This function approximates this term as
     * $ - \nabla \mathbf{u} = \frac{1}{\rho} \frac{\partial rho}{\partial p} \rho \mathbf{g} \cdot \mathbf{u}$
     */
    template <int dim>
    class NewtonStokesIsentropicCompressionTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the right-hand-side term of the Stokes equation
     * that is caused by the variable density in the mass conservation equation.
     * This class approximates this term as
     * $ - \nabla \cdot \mathbf{u} = \frac{1}{\rho} \frac{\partial \rho}{\partial t} + \frac{1}{\rho} \nabla \rho \cdot \mathbf{u}$
     * where the right-hand side velocity is explicitly taken from the last timestep,
     * and the density is taken from a compositional field of the type 'density'.
     */
    template <int dim>
    class NewtonStokesProjectedDensityFieldTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const override;
    };
  }
}

#endif
