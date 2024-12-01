/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/diffusion_dislocation.h>
#include <aspect/utilities.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    DiffusionDislocation<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
        {
          // const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const std::vector<double> composition = in.composition[i];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions(composition,
                                                       this->introspection().chemical_composition_field_indices());

          // densities
          double density = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            {
              // not strictly correct if thermal expansivities are different, since we are interpreting
              // these compositions as volume fractions, but the error introduced should not be too bad.
              const double temperature_factor= (1.0 - thermal_expansivities[j] * (temperature - reference_T));
              density += volume_fractions[j] * densities[j] * temperature_factor;
            }

          // thermal expansivities
          double thermal_expansivity = 0.0;
          for (unsigned int j=0; j < volume_fractions.size(); ++j)
            thermal_expansivity += volume_fractions[j] * thermal_expansivities[j];

          // calculate effective viscosity
          if (in.requests_property(MaterialProperties::viscosity))
            {
              Assert(std::isfinite(in.strain_rate[i].norm()),
                     ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                "not filled by the caller."));
              out.viscosities[i] = diffusion_dislocation.compute_viscosity(pressure, temperature, volume_fractions, in.strain_rate[i]);
            }

          out.densities[i] = density;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = heat_capacity;
          // Thermal conductivity at the given positions. If the temperature equation uses
          // the reference density profile formulation, use the reference density to
          // calculate thermal conductivity. Otherwise, use the real density. If the adiabatic
          // conditions are not yet initialized, the real density will still be used.
          if (this->get_parameters().formulation_temperature_equation ==
              Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile &&
              this->get_adiabatic_conditions().is_initialized())
            out.thermal_conductivities[i] = thermal_diffusivity * heat_capacity *
                                            this->get_adiabatic_conditions().density(in.position[i]);
          else
            out.thermal_conductivities[i] = thermal_diffusivity * heat_capacity * density;
          // Compressibility at the given positions.
          // The compressibility is given as
          // $\frac 1\rho \frac{\partial\rho}{\partial p}$.
          out.compressibilities[i] = 0.0;
          // Pressure derivative of entropy at the given positions.
          out.entropy_derivative_pressure[i] = 0.0;
          // Temperature derivative of entropy at the given positions.
          out.entropy_derivative_temperature[i] = 0.0;
          // Change in composition due to chemical reactions at the
          // given positions. The term reaction_terms[i][c] is the
          // change in compositional field c at point i.
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }

    template <int dim>
    bool
    DiffusionDislocation<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    DiffusionDislocation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Diffusion dislocation");
        {
          Rheology::DiffusionDislocation<dim>::declare_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DiffusionDislocation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Diffusion dislocation");
        {
          reference_T = prm.get_double("Reference temperature");

          // Retrieve the list of composition names
          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
          // Establish that a background field is required here
          compositional_field_names.insert(compositional_field_names.begin(),"background");

          // Make options file for parsing maps to double arrays
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

          // Establish that a background field is required here
          chemical_field_names.insert(chemical_field_names.begin(),"background");

          Utilities::MapParsing::Options options(chemical_field_names, "");
          options.list_of_allowed_keys = compositional_field_names;
          options.property_name = "Densities";
          densities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          options.property_name = "Thermal expansivities";
          thermal_expansivities = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

          // Phenomenological parameters
          thermal_diffusivity = prm.get_double("Thermal diffusivity");
          heat_capacity = prm.get_double("Heat capacity");

          diffusion_dislocation.initialize_simulator (this->get_simulator());
          diffusion_dislocation.parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate | NonlinearDependence::compositional_fields;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DiffusionDislocation,
                                   "diffusion dislocation",
                                   "An implementation of a viscous rheology including diffusion "
                                   "and dislocation creep. "
                                   "Compositional fields can each be assigned individual "
                                   "activation energies, reference densities, thermal expansivities, "
                                   "and stress exponents. The effective viscosity is defined as "
                                   "\n\n"
                                   "$\\eta_{\\text{eff}} = \\left(\\frac{1}{\\eta_{\\text{eff}}^\\text{diff}}+ "
                                   "\\frac{1}{\\eta_{\\text{eff}}^\\text{dis}}\\right)^{-1}$ "
                                   "where "
                                   "$\\eta_{\\text{i}} = \\frac{1}{2} A^{-\\frac{1}{n_i}} d^\\frac{m_i}{n_i} "
                                   "\\dot{\\varepsilon_i}^{\\frac{1-n_i}{n_i}} "
                                   "\\exp\\left(\\frac{E_i^\\ast + PV_i^\\ast}{n_iRT}\\right)$ "
                                   "\n\n"
                                   "where $d$ is grain size, $i$ corresponds to diffusion or dislocation creep, "
                                   "$\\dot{\\varepsilon}$ is the square root of the second invariant of the "
                                   "strain rate tensor, $R$ is the gas constant, $T$ is temperature, "
                                   "and $P$ is pressure. "
                                   "$A_i$ are prefactors, $n_i$ and $m_i$ are stress and grain size exponents "
                                   "$E_i$ are the activation energies and $V_i$ are the activation volumes. "
                                   "\n\n"
                                   "This form of the viscosity equation is commonly used in geodynamic simulations "
                                   "See, for example, Billen and Hirth (2007), G3, 8, Q08012. Significantly, "
                                   "other studies may use slightly different forms of the viscosity equation "
                                   "leading to variations in how specific terms are defined or combined. For "
                                   "example, the grain size exponent should always be positive in the diffusion "
                                   "viscosity equation used here, while other studies place the grain size term "
                                   "in the denominator and invert the sign of the grain size exponent. When "
                                   "examining previous work, one should carefully check how the viscous "
                                   "prefactor and grain size terms are defined. "
                                   " \n\n"
                                   "The ratio of diffusion to dislocation strain rate is found by Newton's "
                                   "method, iterating to find the stress which satisfies the above equations. "
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the parameter file in subsection "
                                   "'Material model/DiffusionDislocation'.")
  }
}
