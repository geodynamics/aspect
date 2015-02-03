/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#ifndef __aspect__model_morency_doin_h
#define __aspect__model_morency_doin_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model based on the rheology described in (Morency and Doin,
     * 2004): Brittle-ductile rheology with a viscosity strongly depending on
     * temperature and composition. Using pseudo-brittle rheology limits the
     * strength of the lithosphere at low temperature.
     *
     * The effective viscosity is defined as the harmonic mean of two stress-
     * dependent viscosity functions: a simple temperature/pressure-dependent,
     * non-newtonian viscosity, and a more strongly stress-dependent,
     * "plastic" viscosity.
     *
     * @f[ v_{eff}^v = B \left(\frac{\dot{\varepsilon}}{\dot{\varepsilon}_{ref}}\right)^{-1+1/n_v}
     * exp\left(\frac{E_a +V_a \rho_m g z}{n_v R T}\right) @f]
     * @f[ v_{eff}^p = (\tau_0 + \gamma \rho_m g z) \left( \frac{\dot{\varepsilon}^{-1+1/n_p}}
     * {\dot{\varepsilon}_{ref}^{1/n_p}} \right) @f]
     * @f[ v_{eff} = \left(\frac{1}{v_{eff}^v}+\frac{1}{v_{eff}^p}\right)^{-1} @f]
     *
     * Where $v_{eff}$ is the effective viscosity, $B$ is a scaling constant,
     * $\dot{\varepsilon}$ is related to the second invariant of the strain
     * rate tensor, $\dot{\varepsilon}_{ref}$ is a reference strain rate,
     * $n_v$ and $n_p$ are stress exponents, $E_a$ is the activation energy,
     * $V_a$ is the activation volume, $\rho_m$ is the mantle density, $R$ is
     * the gas constant, $T$ is temperature, $\tau_0$ is the cohestive
     * strength of rocks at the surface, $\gamma$ is a coefficient of yield
     * stress increase with depth, and $z$ is depth.
     *
     * Several model parameters (reference densities, activation energies,
     * thermal expansivities, and stress exponents-- both viscous and plastic)
     * can be defined per-compositional field. If a list of values is given
     * for any of these parameters, the weighted sum of the values based on
     * volume fractions of the compositional fields is used in their place. If
     * only one value is given for any of these parameters, all compositions
     * are assigned the same value. The first value in the list is the value
     * assigned to "background mantle" (regions where the sum of the
     * compositional fields is < 1.0).
     *
     * For more on the material model and its applications, see: Morency, C.,
     * and M‐P. Doin. "Numerical simulations of the mantle lithosphere
     * delamination." Journal of Geophysical Research: Solid Earth
     * (1978–2012) 109.B3 (2004)
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MorencyDoin : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        typedef typename aspect::MaterialModel::Interface<dim>::MaterialModelInputs MaterialModelInputs;
        typedef typename aspect::MaterialModel::Interface<dim>::MaterialModelOutputs MaterialModelOutputs;

        virtual void evaluate(const MaterialModelInputs &in, MaterialModelOutputs &out) const;


        /**
         * Return true if the viscosity() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
         *
         * This function must return false for all possible arguments if the
         * is_compressible() function returns false.
         */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the specific_heat() function returns something that
         * may depend on the variable identifies by the argument.
         */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
         */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
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
        /**
         * From multicomponent material model: From a list of compositional
         * fields of length N, we come up with an N+1 length list that which
         * also includes the fraction of ``background mantle''. This list
         * should sum to one, and is interpreted as volume fractions.  If the
         * sum of the compositional_fields is greater than one, we assume that
         * there is no background mantle (i.e., that field value is zero).
         * Otherwise, the difference between the sum of the compositional
         * fields and 1.0 is assumed to be the amount of background mantle.
         */
        std::vector<double> compute_volume_fractions(
          const std::vector<double> &compositional_fields) const;
        std::vector<double> densities;
        std::vector<double> activation_energies;
        std::vector<double> thermal_expansivities;
        std::vector<double> nvs; //Stress exponent, viscous rheology
        std::vector<double> nps;//Stress exponent, plastic rheology
        double thermal_diffusivity;
        double gamma; // Coefficient of yield stress increase with depth
        double heat_capacity;
        double activation_volume;
        double ref_strain_rate;
        double B; //Preexponential constant in the viscous rheology law B
        double tau_0; //cohesive strength of rocks at the surface
        double reference_T;
        double min_strain_rate;

        double min_visc;
        double max_visc;
        double veff_coefficient;
        double ref_visc;
    };

  }
}

#endif
