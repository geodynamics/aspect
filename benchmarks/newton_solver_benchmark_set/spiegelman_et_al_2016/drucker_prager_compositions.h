/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#ifndef __aspect__model_drucker_prager_compositions_h
#define __aspect__model_drucker_prager_compositions_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/newton.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * The same material model as Drucker Prager, but this one supports multiple
     * compositions. Designed to run the Spiegelman et al. 2016 benchmark.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * The viscosity is computed according to the Drucker Prager frictional
     * plasticity criterion based on a user-defined internal angle of friction $\phi$
     * and cohesion $C$. In 3D:
     * $\sigma_y = \frac{6 C \cos(\phi)}{\sqrt(3) (3+\sin(\phi))} +
     * \frac{6 P \sin(\phi)}{\sqrt(3) (3+\sin(\phi))}$,
     * where $P$ is the pressure.
     * See for example Zienkiewicz, O. C., Humpheson, C. and Lewis, R. W. (1975),
     * Géotechnique 25, No. 4, 671-689.
     * With this formulation we circumscribe instead of inscribe the Mohr Coulomb
     * yield surface.
     * In 2D the Drucker Prager yield surface is the same
     * as the Mohr Coulomb surface:
     * $\sigma_y = P \sin(\phi) + C \cos(\phi)$.
     * Note that in 2D for $\phi=0$, these criteria
     * revert to the von Mises criterion (no pressure dependence).
     * See for example Thieulot, C. (2011), PEPI 188, 47-68.
     *
     * Note that we enforce the pressure to be positive in the computation of
     * the yield strength by replacing it with
     * a zero value whenever it is negative to prevent negative
     * yield strengths and viscosities.
     * We then use the computed yield strength to scale back the viscosity on
     * to the yield surface using the Viscosity Rescaling Method described in
     * Kachanov, L. M. (2004), Fundamentals of the Theory of Plasticity,
     * Dover Publications, Inc.
     *
     * To avoid numerically unfavourably large (or even negative) viscosity ranges,
     * we cut off the viscosity with a user-defined minimum and maximum viscosity:
     * $\eta_eff = \frac{1}{\frac{1}{\eta_min + \eta}+\\
     * \frac{1}{\eta_max}}$.
     *
     * Note that this model uses the formulation that assumes an incompressible
     * medium despite the fact that the density follows the law
     * $\rho(T)=\rho_0(1-\beta(T-T_{\text{ref}}))$.
     *
     *
     * @ingroup MaterialModels
     */

    template <int dim>
    class DruckerPragerCompositions : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        std::vector<double> compute_volume_fractions( const std::vector<double> &compositional_fields) const;

        double compute_second_invariant(const SymmetricTensor<2,dim> strain_rate, const double min_strain_rate) const;

        double compute_viscosity(const double edot_ii,const double pressure,const int comp, const double prefactor,const bool regularize, const double min_visc, const double max_visc) const;

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
        *
        * This material model is incompressible.
         */
        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        double reference_T;
        std::vector<double> reference_T_list;

        /**
         * Defining a minimum strain rate stabilizes the viscosity calculation,
         * which involves a division by the strain rate. Units: $1/s$.
         */
        bool use_deviator_of_strain_rate;
        std::vector<double> min_strain_rate;
        std::vector<double> min_visc;
        std::vector<double> max_visc;
        std::vector<double> veff_coefficient;
        double ref_visc;
        double reference_compressibility;
        std::vector<double> ref_visc_list;

        std::vector<double> thermal_diffusivity;
        std::vector<double> heat_capacity;


        std::vector<double> densities;
        std::vector<double> thermal_expansivities;

        std::vector<double> cohesion;
        std::vector<double> phi;

        std::vector<double> prefactor;
        std::vector<double> stress_exponent;

        unsigned int n_fields;
        // averaging parameters
        double viscosity_averaging_p;

        bool use_analytical_derivative;

        //optimize variables
        const double sqrt_3 = std::sqrt(3);
        const double sqrt_half = std::sqrt(0.5);
        std::vector<double> sin_phi;
        std::vector<double> cos_phi;



    };

  }
}

#endif
