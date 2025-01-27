/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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

#include <aspect/material_model/thermal_conductivity/tosi_stackhouse.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
      template <int dim>
      void
      TosiStackhouse<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                                     MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        const unsigned int n_evaluation_points = in.n_evaluation_points();
        for (unsigned int i=0; i<n_evaluation_points; ++i)
          {
            // Find the conductivity layer that corresponds to the depth of the evaluation point.
            const double depth = this->get_geometry_model().depth(in.position[i]);
            const unsigned int layer_index = std::distance(conductivity_transition_depths.begin(),
                                                           std::lower_bound(conductivity_transition_depths.begin(),conductivity_transition_depths.end(), depth));

            const double p_dependence = reference_thermal_conductivities[layer_index] + conductivity_pressure_dependencies[layer_index] * in.pressure[i];

            // Make reasonably sure we will not compute any invalid values due to the temperature-dependence.
            // Since both the temperature-dependence and the saturation term scale with (Tref/T), we have to
            // make sure we can compute the square of this number. If the temperature is small enough to
            // be close to yielding NaN values, the conductivity will be set to the maximum value anyway.
            const double T = std::max(in.temperature[i], std::sqrt(std::numeric_limits<double>::min()) * conductivity_reference_temperatures[layer_index]);
            const double T_dependence = std::pow(conductivity_reference_temperatures[layer_index] / T, conductivity_exponents[layer_index]);

            // Function based on the theory of Roufosse and Klemens (1974) that accounts for saturation.
            // For the Tosi formulation, the scaling should be zero so that this term is 1.
            double saturation_function = 1.0;
            if (1./T_dependence > 1.)
              saturation_function = (1. - saturation_scaling[layer_index])
                                    + saturation_scaling[layer_index] * (2./3. * std::sqrt(T_dependence) + 1./3. * 1./T_dependence);

            out.thermal_conductivities[i] = std::min(p_dependence * saturation_function * T_dependence, maximum_conductivity);
          }
      }



      template <int dim>
      void
      TosiStackhouse<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Thermal conductivity transition depths", "410000, 520000, 660000",
                           Patterns::List(Patterns::Double (0.)),
                           "A list of depth values that indicate where the transitions between "
                           "the different conductivity parameter sets should occur (in most cases, "
                           "these will be the depths of major phase transitions). "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Reference thermal conductivities", "2.47, 3.81, 3.52, 4.9",
                           Patterns::List(Patterns::Double (0.)),
                           "A list of base values of the thermal conductivity for each of the "
                           "horizontal layers. Pressure- and temperature-dependence will be applied "
                           "on top of this base value, according to the parameters 'Pressure "
                           "dependencies of thermal conductivity' and 'Reference temperatures "
                           "for thermal conductivity'. "
                           "Units: \\si{\\watt\\per\\meter\\per\\kelvin}");
        prm.declare_entry ("Pressure dependencies of thermal conductivity", "3.3e-10, 3.4e-10, 3.6e-10, 1.05e-10",
                           Patterns::List(Patterns::Double ()),
                           "A list of values that determine the linear scaling of the "
                           "thermal conductivity with pressure. "
                           "Units: \\si{\\watt\\per\\meter\\per\\kelvin\\per\\pascal}.");
        prm.declare_entry ("Reference temperatures for thermal conductivity", "300, 300, 300, 1200",
                           Patterns::List(Patterns::Double (0.)),
                           "A list of values of reference temperatures used to determine "
                           "the temperature-dependence of the thermal conductivity. "
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Thermal conductivity exponents", "0.48, 0.56, 0.61, 1.0",
                           Patterns::List(Patterns::Double (0.)),
                           "A list of exponents in the temperature-dependent term of the "
                           "conductivity formulation. Note that this "
                           "exponent is not used (and should have a value of 1) in the "
                           "formulation of Stackhouse et al. (2015). "
                           "Units: none.");
        prm.declare_entry ("Saturation prefactors", "0, 0, 0, 1",
                           Patterns::List(Patterns::Double (0., 1.)),
                           "A list of values that indicate how a given layer "
                           "should take into account the effects "
                           "of saturation on the temperature-dependence of the thermal "
                           "conductivity. This factor is multiplied with a saturation function "
                           "based on the theory of Roufosse and Klemens, 1974. A value of 1 "
                           "reproduces the formulation of Stackhouse et al. (2015), a value of "
                           "0 reproduces the formulation of Tosi et al., (2013). "
                           "Units: none.");
        prm.declare_entry ("Maximum thermal conductivity", "1000",
                           Patterns::Double (0.),
                           "The maximum thermal conductivity that is allowed in the "
                           "model. Larger values will be cut off.");
      }



      template <int dim>
      void
      TosiStackhouse<dim>::parse_parameters (ParameterHandler &prm)
      {
        conductivity_transition_depths = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("Thermal conductivity transition depths")));
        const unsigned int n_conductivity_layers = conductivity_transition_depths.size() + 1;

        AssertThrow (std::is_sorted(conductivity_transition_depths.begin(), conductivity_transition_depths.end()),
                     ExcMessage("The list of 'Thermal conductivity transition depths' must "
                                "be sorted such that the values increase monotonically."));

        reference_thermal_conductivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference thermal conductivities"))),
                                                                                   n_conductivity_layers,
                                                                                   "Reference thermal conductivities");
        conductivity_pressure_dependencies = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Pressure dependencies of thermal conductivity"))),
                                                                                     n_conductivity_layers,
                                                                                     "Pressure dependencies of thermal conductivity");
        conductivity_reference_temperatures = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Reference temperatures for thermal conductivity"))),
                                                                                      n_conductivity_layers,
                                                                                      "Reference temperatures for thermal conductivity");
        conductivity_exponents = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal conductivity exponents"))),
                                                                         n_conductivity_layers,
                                                                         "Thermal conductivity exponents");
        saturation_scaling = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Saturation prefactors"))),
                                                                     n_conductivity_layers,
                                                                     "Saturation prefactors");
        maximum_conductivity = prm.get_double ("Maximum thermal conductivity");
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace ThermalConductivity
    {
#define INSTANTIATE(dim) \
  template class TosiStackhouse<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
