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


#include <aspect/material_model/nondimensional.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Nondimensional<dim>::initialize ()
    {
      Assert(this->get_parameters().formulation_temperature_equation
             == Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile,
             ExcMessage("The Nondimensional material model can only work with a "
                        "temperature formulation that is based on the reference "
                        "profile."));
      Assert((this->get_parameters().formulation_mass_conservation
              == Parameters<dim>::Formulation::MassConservation::incompressible
              && !compressible)
             ||
             (this->get_parameters().formulation_mass_conservation
              == Parameters<dim>::Formulation::MassConservation::reference_density_profile
              && compressible)
             ||
             (this->get_parameters().formulation_mass_conservation
              == Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile
              && compressible)
             ,
             ExcMessage("The Nondimensional material model can only work with a "
                        "mass formulation that is incompressible or based on "
                        "the reference profile and the Di has to be chosen "
                        "accordingly."));
    }


    template <int dim>
    void
    Nondimensional<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature_deviation = in.temperature[i] - this->get_adiabatic_conditions().temperature(position);
          const double pressure_deviation = in.pressure[i] - this->get_adiabatic_conditions().pressure(position);

          const double depth = this->get_geometry_model().depth(position);

          out.viscosities[i] = (compressible ? Di : 1.0) / Ra
                               * std::exp( -exponential_viscosity_temperature_prefactor * temperature_deviation
                                           + exponential_viscosity_depth_prefactor * depth);

          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = 1.0;
          out.thermal_expansion_coefficients[i] = compressible ? Di : 1.0;

          out.densities[i] = reference_rho
                             * std::exp(depth * Di/gamma)
                             * (1.0
                                - out.thermal_expansion_coefficients[i] * temperature_deviation
                                + (tala ? 0.0 : 1.0) * Di * gamma * pressure_deviation);

          out.compressibilities[i] = 0.0;
          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }


    template <int dim>
    double
    Nondimensional<dim>::
    reference_viscosity () const
    {
      return compressible ? (Di/Ra) : (1.0/Ra);
    }



    template <int dim>
    bool
    Nondimensional<dim>::
    is_compressible () const
    {
      return compressible;
    }



    template <int dim>
    void
    Nondimensional<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Nondimensional model");
        {
          prm.declare_entry ("Reference density", "1.0",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Ra", "1e4",
                             Patterns::Double (0),
                             "Rayleigh number Ra");
          prm.declare_entry ("Di", "0.0",
                             Patterns::Double (0),
                             "Dissipation number. Pick 0.0 for incompressible "
                             "computations.");
          prm.declare_entry ("gamma", "1.0",
                             Patterns::Double (0),
                             "Grueneisen parameter");
          prm.declare_entry ("Reference specific heat", "1.0",
                             Patterns::Double (0),
                             "The value of the specific heat $C_p$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Viscosity temperature prefactor", "0.0",
                             Patterns::Double (0),
                             "Exponential temperature prefactor for viscosity.");
          prm.declare_entry ("Viscosity depth prefactor", "0.0",
                             Patterns::Double (0),
                             "Exponential depth prefactor for viscosity.");
          prm.declare_entry ("Use TALA", "false",
                             Patterns::Bool (),
                             "Whether to use the TALA instead of the ALA "
                             "approximation.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Nondimensional<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Nondimensional model");
        {
          reference_rho   = prm.get_double ("Reference density");
          exponential_viscosity_temperature_prefactor = prm.get_double ("Viscosity temperature prefactor");
          exponential_viscosity_depth_prefactor = prm.get_double ("Viscosity depth prefactor");
          Di              = prm.get_double ("Di");
          Ra              = prm.get_double ("Ra");
          gamma           = prm.get_double ("gamma");
          tala            = prm.get_bool ("Use TALA");
          reference_specific_heat = prm.get_double ("Reference specific heat");

          compressible = (Di!=0.0);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature;
      this->model_dependence.density = NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Nondimensional,
                                   "nondimensional",
                                   "A material model for nondimensionalized "
                                   "computations for compressible or incompressible "
                                   "computations defined through Rayleigh number \text{Ra} "
                                   "and Dissipation number Di. This model is made "
                                   "to be used with the Boussinesq, ALA, or TALA "
                                   "formulation."
                                   "\n\n"
                                   "The viscosity is defined as \\[\\eta = \\text{Di} / \\text{Ra} \\cdot \\exp(-b T' + c z)\\] "
                                   "where $T'$ is the temperature variation from the "
                                   "adiabatic temperature, $z$ is the depth, "
                                   "$b$ is given by ``Viscosity temperature prefactor'', "
                                   "and $c$ by ``Viscosity depth prefactor''. If "
                                   "$\\text{Di}$ is zero, it will be replaced by 1.0 in $\\eta$."
                                   "\n\n"
                                   "The density is defined as \\[\\rho = \\exp(\\text{Di}/\\gamma \\cdot z) "
                                   " (1.0 - \\alpha T' + \\text{Di} \\gamma p'),\\] where "
                                   "$\\alpha=\text{Di}$ is the thermal expansion coefficient, "
                                   "$\\gamma$ is the Grueneisen parameter, and $p'$ is "
                                   "the pressure variation from the adiabatic "
                                   "pressure. The pressure dependent term is not present "
                                   "if ``TALA'' is enabled.")
  }
}
