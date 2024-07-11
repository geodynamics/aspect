/*
  Copyright (C) 2017 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_simulator_assemblers_stokes_h
#define _aspect_simulator_assemblers_stokes_h


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Assemblers
  {
    /**
     * A class containing the functions to assemble the Stokes preconditioner.
     */
    template <int dim>
    class StokesPreconditioner : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * A class containing the functions to assemble the compressible adjustment
     * to the Stokes preconditioner.
     */
    template <int dim>
    class StokesCompressiblePreconditioner : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Stokes equation for the current cell.
     */
    template <int dim>
    class StokesIncompressibleTerms : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;

        /**
         * Create AdditionalMaterialOutputsStokesRHS if we need to do so.
         */
        void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * This class assembles the term that arises in the viscosity term of Stokes matrix for
     * compressible models, because the divergence of the velocity is not longer zero.
     */
    template <int dim>
    class StokesCompressibleStrainRateViscosityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the compressibility term of the Stokes equation
     * that is caused by the compressibility in the mass conservation equation.
     * It uses an approximation that involves the reference density profile, and
     * includes this term explicitly in the right-hand side vector to preserve
     * the symmetry of the matrix.
     * This class approximates this term as
     * $- \nabla \cdot \mathbf{u} = \frac{1}{\rho^{\ast}} \frac{\partial rho}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u}$
     */
    template <int dim>
    class StokesReferenceDensityCompressibilityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the compressibility term of the Stokes equation
     * that is caused by the compressibility in the mass conservation equation.
     * It uses an approximation that involves the reference density profile, and
     * includes this term implicitly in the matrix,
     * which is therefore not longer symmetric.
     * This class approximates this term as
     * $ - \nabla \cdot \mathbf{u} - \frac{1}{\rho^{\ast}} \frac{\partial rho{^\ast}}{\partial z} \frac{\mathbf{g}}{||\mathbf{g}||} \cdot \mathbf{u} = 0$
     */
    template <int dim>
    class StokesImplicitReferenceDensityCompressibilityTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the right-hand-side term of the Stokes equation
     * that is caused by the compressibility in the mass conservation equation.
     * This class approximates this term as
     * $ - \nabla \cdot \mathbf{u} = \kappa \rho \mathbf{g} \cdot \mathbf{u}$
     * where $\kappa$ is the compressibility provided by the material model,
     * which is frequently computed as
     * $\kappa = \frac{1}{\rho} \frac{\partial rho}{\partial p}$.
     */
    template <int dim>
    class StokesIsentropicCompressionTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
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
    class StokesProjectedDensityFieldTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const override;
    };


    /**
     * This class assembles the right-hand-side term of the Stokes equation
     * that is caused by the compression based on the changes in (hydrostatic)
     * pressure and temperature in the mass conservation equation.
     *
     * This class approximates this term as
     * $ -\nabla \cdot \mathbf{u} = \left( \kappa \rho \textbf{g} - \alpha \nabla T \right) \cdot \textbf{u}$
     *
     * where $\frac{1}{\rho} \frac{\partial \rho}{\partial p} = \kappa$ is the compressibility,
     * $- \frac{1}{\rho}\frac{\partial \rho}{\partial T} = \alpha$ is the thermal expansion coefficient,
     * and both are defined in the material model.
     */
    template <int dim>
    class StokesHydrostaticCompressionTerm : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class assembles the right-hand-side terms that are used to weakly
     * prescribe the boundary tractions.
     */
    template <int dim>
    class StokesBoundaryTraction : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };

    /**
     * This class computes the local pressure shape function integrals that
     * are later used to make the Stokes equations compatible to its right hand
     * side. For more information why this is necessary see Section 3.2.2 of
     * Heister et al. (2017), "High Accuracy Mantle Convection Simulation
     * through Modern Numerical Methods. II: Realistic Models and Problems."
     */
    template <int dim>
    class StokesPressureRHSCompatibilityModification : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
                internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const override;
    };
  }
}


#endif
