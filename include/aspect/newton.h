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


#ifndef _aspect__newton_h
#define _aspect__newton_h

#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/material_model/interface.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  using namespace dealii;

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
         * The derivatives of the viscosities
         */
        std::vector<double> viscosity_derivative_wrt_pressure;
        std::vector<SymmetricTensor<2,dim> > viscosity_derivative_wrt_strain_rate;

    };
  }

  /**
   * A Class which can declare and parse parameters and creates
   * material model outputs for the Newton solver.
   */
  template <int dim>
  class NewtonHandler: public SimulatorAccess<dim>
  {
    public:
      /**
       * Create an additional material model output object that contains
       * the additional output variables (the derivatives) needed for the
       * Newton solver.
       */
      static void create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output);

      /**
       * Return the Newton derivative scaling factor used for scaling the
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
      double get_newton_derivative_scaling_factor() const;

      /**
       * Set the Newton derivative scaling factor used for scaling the
       * derivative part of the Newton Stokes solver in the assembly.
       *
       * See the get_newton_derivative_scaling_factor() function for an
       * explanation of the purpose of this factor.
       */
      void set_newton_derivative_scaling_factor(const double newton_derivative_scaling_factor);


    private:
      /**
       * A scaling factor for those terms of the Newton matrix that
       * result from the linearization of the viscosity.
       *
       * See the get_newton_derivative_scaling_factor() function for an
       * explanation of the purpose of this factor.
       */
      double newton_derivative_scaling_factor;
  };

  namespace Assemblers
  {
    /**
      * A base class for the definition of assemblers that implement the
      * linear system terms for the NewtonStokes solver scheme.
      */
    template <int dim>
    class NewtonInterface : public aspect::Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual ~NewtonInterface () {};

        /**
         * Attach Newton outputs. Since most Newton assemblers require the
         * material model derivatives they are created in this base class
         * already.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };

    /**
     * This class assembles the terms of the Newton Stokes preconditioner matrix for the current cell.
     */
    template <int dim>
    class NewtonStokesPreconditioner : public NewtonInterface<dim>
    {
      public:
        virtual ~NewtonStokesPreconditioner () {};

        void
        execute (internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Newton Stokes system for the current cell.
     */
    template <int dim>
    class NewtonStokesIncompressibleTerms : public NewtonInterface<dim>
    {
      public:
        virtual ~NewtonStokesIncompressibleTerms () {};

        void
        execute (internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
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
        virtual ~NewtonStokesCompressibleStrainRateViscosityTerm () {};

        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
    };

    /**
     * This class assembles the right-hand-side term of the Newton Stokes system
     * that is caused by the compressibility in the mass conservation equation.
     * This function approximates this term as
     * $- \nabla \mathbf{u} = \frac{1}{\rho} * \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u}$
     */
    template <int dim>
    class NewtonStokesReferenceDensityCompressibilityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual ~NewtonStokesReferenceDensityCompressibilityTerm () {};

        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
    };

    /**
     * This class assembles the compressibility term of the Newton Stokes system
     * that is caused by the compressibility in the mass conservation equation.
     * It includes this term implicitly in the matrix,
     * which is therefore not longer symmetric.
     * This function approximates this term as
     * $ - \nabla \mathbf{u} - \frac{1}{\rho} * \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u} = 0$
     */
    template <int dim>
    class NewtonStokesImplicitReferenceDensityCompressibilityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual ~NewtonStokesImplicitReferenceDensityCompressibilityTerm () {};

        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
    };

    /**
     * This class assembles the right-hand-side term of the Newton Stokes system
     * that is caused by the compressibility in the mass conservation equation.
     * This function approximates this term as
     * $ - \nabla \mathbf{u} = \frac{1}{\rho} * \frac{\partial rho}{\partial p} \rho \mathbf{g} \cdot \mathbf{u}$
     */
    template <int dim>
    class NewtonStokesIsothermalCompressionTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual ~NewtonStokesIsothermalCompressionTerm () {};

        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>  &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
    };
  }
}

#endif
