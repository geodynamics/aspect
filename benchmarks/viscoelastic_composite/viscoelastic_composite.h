/*
  Copyright (C) 2020 - 2026 by the authors of the ASPECT code.

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
#ifndef _aspect_material_model_viscoelastic_composite_h
#define _aspect_material_model_viscoelastic_composite_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/material_model/equation_of_state/linearized_incompressible.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model implementing a generalized Maxwell (Wiechert) rheology:
     * an arbitrary number of Maxwell arms (spring + dashpot in series) connected
     * in parallel. Since all arms share the same total strain rate, their
     * stress contributions -- and effective viscosities -- add linearly; see
     * evaluate() and fill_elastic_outputs() in the .cc file.
     *
     * Loosely motivated by the two-Maxwell-arm crustal rheology used in
     * Head et al. (2019), "The Influence of Viscoelastic Crustal Rheologies
     * on Volcanic Ground Deformation: Insights From Models of Pressure and
     * Volume Change".
     */
    template <int dim>
    class ViscoelasticComposite : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * This material model is incompressible.
         */
        bool is_compressible () const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        void
        create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:
        /**
         * Compute the total (summed) elastic force term across all Maxwell
         * arms and attach it to the ElasticOutputs additional output object
         * on `out`, so it reaches the Stokes assembly.
         */
        void fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                   const std::vector<double> &average_elastic_shear_moduli,
                                   MaterialModel::MaterialModelOutputs<dim> &out) const;

        double minimum_strain_rate;

        double elastic_timestep () const;

        /**
         * The elasticity rheology module shared across all arms. Currently
         * assumes a single elastic timestep and rotation scheme is applied
         * uniformly; per-arm elastic behavior (relaxation time) comes from
         * pairing elastic_shear_moduli[n] with viscosities[n], not from
         * separate Elasticity objects per arm.
         */
        Rheology::Elasticity<dim> elastic_rheology;

        /**
         * Viscosity of each Maxwell arm's dashpot, read from the parameter
         * file. One entry per arm, ordered to match elastic_shear_moduli.
         */
        std::vector<double> viscosities;

        /**
         * Elastic shear modulus of each Maxwell arm's spring, read from the
         * parameter file. One entry per arm, ordered to match viscosities.
         */
        std::vector<double> elastic_shear_moduli;

        /**
         * Number of parallel Maxwell arms in the composite model. Read from
         * the parameter file in parse_parameters(); must equal the number
         * of entries in `viscosities` and `elastic_shear_moduli` after
         * possibly_extend_from_1_to_N is applied.
         */
        unsigned int number_of_maxwell_arms;

        /**
         * Number of independent components of a symmetric rank-2 tensor in
         * `dim` dimensions (3 in 2D, 6 in 3D). Used in fill_elastic_outputs()
         * to compute compositional field offsets for each arm's stored
         * stress history (advected + old), which are laid out per arm as
         * 2 * n_independent_components contiguous fields.
         */
        static const unsigned int n_independent_components = SymmetricTensor<2,dim>::n_independent_components;

        EquationOfState::LinearizedIncompressible<dim> equation_of_state;
    };
  }
}

#endif
