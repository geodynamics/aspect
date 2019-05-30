/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class FrankKamenetskii : public MaterialModel::Simple<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
        virtual void parse_parameters(ParameterHandler &prm);
        static void declare_parameters(ParameterHandler &prm);
        virtual double reference_viscosity () const;

      private:
        double viscosity_reference_T;
        double viscosity_variation;
        double eta;
        double composition_viscosity_prefactor;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    FrankKamenetskii<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // First, we use the material descriptions of the 'simple' material model to fill all of the material
      // model outputs. We will then overwrite selected properties.
      Simple<dim>::evaluate(in, out);

      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          // Calculate the temperature perturbation for the viscosity law.
          const double delta_temp = in.temperature[i]-viscosity_reference_T;

          // Calculate the activation energy, E, from viscosity variation, delta_eta.
          // This follows the method from Zhong et al., 2008.
          const double activation_energy = std::log(viscosity_variation);

          // Calculate the temperature dependence for the viscosity law.
          const double temperature_dependence = std::exp(-activation_energy * delta_temp);

          // Calculate the temperature-dependent viscosity.
          out.viscosities[i] = ((composition_viscosity_prefactor != 1.0) && (in.composition[i].size()>0))
                               ?
                               // Geometric interpolation
                               std::pow(10.0, ((1-in.composition[i][0]) * std::log10(eta * temperature_dependence)
                                               + in.composition[i][0] * std::log10(eta * composition_viscosity_prefactor * temperature_dependence)))
                               :
                               temperature_dependence * eta;

        }
    }

    template <int dim>
    double
    FrankKamenetskii<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    void
    FrankKamenetskii<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the constant viscosity $\\eta_0$. This viscosity may be "
                             "modified by both temperature and compositional dependencies. Units: $kg/m/s$.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Viscosity reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature used in the viscosity formula. "
                             "Units $K$.");
          prm.declare_entry ("Viscosity variation", "1.0",
                             Patterns::Double (1.0),
                             "Controls the temperature dependence of viscosity in the "
                             "model. Dimensionless exponent. Note that the minimum "
                             "value of this parameter is 1.0, which is equivalent to "
                             "using a constant viscosity throughout the domain.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    FrankKamenetskii<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple model");
        {
          eta                             = prm.get_double ("Viscosity");
          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          viscosity_reference_T           = prm.get_double ("Viscosity reference temperature");
          viscosity_variation             = prm.get_double ("Viscosity variation");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;

      if (viscosity_variation != 1.0)
        this->model_dependence.viscosity |= NonlinearDependence::temperature;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(FrankKamenetskii,
                                   "frank-kamenetskii",
                                   "A simple material model that is like the "
                                   "'Simple' model, but uses a Frank-Kamenetskii rheology. "
                                   "Specifically, formulae for density and viscosity are taken "
                                   "from Zhong et al., 2008. This material model allows the "
                                   "user to use different reference temperature values for "
                                   "density and viscosity.")
  }
}
