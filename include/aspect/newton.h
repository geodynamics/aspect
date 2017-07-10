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

#include <aspect/assembly.h>
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
  };

  namespace Assemblers
  {
    /**
     * A class containing the functions to assemble the different terms of the Newton Stokes system.
     */
    template <int dim>
    class NewtonStokesAssembler : public aspect::internal::Assembly::Assemblers::AssemblerBase<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * This function assembles the terms of the Newton Stokes preconditioner matrix for the current cell.
         */
        void
        preconditioner (const double                                             pressure_scaling,
                        internal::Assembly::Scratch::StokesPreconditioner<dim>  &scratch,
                        internal::Assembly::CopyData::StokesPreconditioner<dim> &data,
                        const Parameters<dim> &parameters) const;

        /**
         * This function assembles the terms for the matrix and right-hand-side of the incompressible
         * Newton Stokes system for the current cell.
         */
        void
        incompressible_terms (const double                                     pressure_scaling,
                              const bool                                       rebuild_stokes_matrix,
                              internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                              internal::Assembly::CopyData::StokesSystem<dim> &data,
                              const Parameters<dim> &parameters) const;

        /**
         * This function assembles the term that arises in the viscosity term of Newton Stokes matrix for
         * compressible models, because the divergence of the velocity is no longer zero.
         */
        void
        compressible_strain_rate_viscosity_term (const double                                     pressure_scaling,
                                                 const bool                                       rebuild_stokes_matrix,
                                                 internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                 internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                 const Parameters<dim> &parameters) const;

        /**
         * This function assembles the right-hand-side term of the Newton Stokes system
         * that is caused by the compressibility in the mass conservation equation.
         * This function approximates this term as
         * $- \nabla \mathbf{u} = \frac{1}{\rho} * \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u}$
         */
        void
        reference_density_compressibility_term (const double                                     pressure_scaling,
                                                const bool                                       rebuild_stokes_matrix,
                                                internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                const Parameters<dim> &parameters) const;

        /**
         * This function assembles the compressibility term of the Newton Stokes system
         * that is caused by the compressibility in the mass conservation equation.
         * It includes this term implicitly in the matrix,
         * which is therefore not longer symmetric.
         * This function approximates this term as
         * $ - \nabla \mathbf{u} - \frac{1}{\rho} * \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u} = 0$
         */
        void
        implicit_reference_density_compressibility_term (const double                                     pressure_scaling,
                                                         const bool                                       rebuild_stokes_matrix,
                                                         internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                                         internal::Assembly::CopyData::StokesSystem<dim> &data,
                                                         const Parameters<dim> &parameters) const;

        /**
         * This function assembles the right-hand-side term of the Newton Stokes system
         * that is caused by the compressibility in the mass conservation equation.
         * This function approximates this term as
         * $ - \nabla \mathbf{u} = \frac{1}{\rho} * \frac{\partial rho}{\partial p} \rho \mathbf{g} \cdot \mathbf{u}$
         */
        void
        isothermal_compression_term (const double                                     pressure_scaling,
                                     const bool                                       rebuild_stokes_matrix,
                                     internal::Assembly::Scratch::StokesSystem<dim>  &scratch,
                                     internal::Assembly::CopyData::StokesSystem<dim> &data,
                                     const Parameters<dim> &parameters) const;

        /**
         * Attach derivatives outputs.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };
  }
}

#endif
