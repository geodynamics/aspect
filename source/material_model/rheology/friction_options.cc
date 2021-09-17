/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#include <aspect/material_model/rheology/friction_options.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      template <int dim>
      double
      FrictionOptions<dim>::
      compute_dependent_friction_angle(const double current_edot_ii,
                                       const unsigned int j,  // j is from a for-loop over volume_fractions.size()
                                       double current_friction) const
      {

        switch (friction_dependence_mechanism)
          {
            case independent:
            {
              break;
            }
            case dynamic_friction:
            {
              // Calculate effective steady-state friction coefficient.
              // This is based on the former material model dynamic friction.
              // The formula below is equivalent to the equation 13 in van Dinther et al., (2013, JGR).
              // Although here the dynamic friction coefficient is directly specified. In addition,
              // we also use a reference strain rate in place of a characteristic
              // velocity divided by local element size. This reference strain rate is called
              // the dynamic characteristic strain rate and is used to see what value between
              // dynamic and static angle of internal friction should be used.
              // Furthermore a smoothness exponent X is added, which determines whether the
              // friction vs strain rate curve is rather step-like or more gradual.
              // mu  = mu_d + (mu_s - mu_d) / ( (1 + strain_rate_dev_inv2/dynamic_characteristic_strain_rate)^X );
              // Angles of friction are used in radians within ASPECT. The coefficient
              // of friction (mu) is the tangent of the internal angle of friction, hence convergence is needed.
              const double mu = (std::tan(dynamic_angles_of_internal_friction[j])
                                 + (std::tan(current_friction) - std::tan(dynamic_angles_of_internal_friction[j]))
                                 / (1. + std::pow((current_edot_ii / dynamic_characteristic_strain_rate),
                                                  dynamic_friction_smoothness_exponent)));
              current_friction = std::atan (mu);
              Assert((mu < 1) && (0 < current_friction) && (current_friction <= 1.6), ExcMessage(
                       "Something is wrong with the tan/atan conversion of friction coefficient to friction angle in RAD."));
              break;
            }
          }
        return current_friction;
      }



      template <int dim>
      FrictionDependenceMechanism
      FrictionOptions<dim>::
      get_friction_dependence_mechanism() const
      {
        return friction_dependence_mechanism;
      }



      template <int dim>
      void
      FrictionOptions<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Friction dependence mechanism", "none",
                           Patterns::Selection("none|dynamic friction"),
                           "Whether to make the friction angle dependent of strain rate. This rheology "
                           "is intended to be used together with the visco-plastic rheology model."
                           "\n\n"
                           "\\item ``none'': No dependence of the friction angle is applied. "
                           "\n\n"
                           "\\item ``dynamic friction'': The friction angle is rate dependent."
                           "When dynamic angles of friction are specified, "
                           "the friction angle will be weakened for high strain rates with: "
                           "$\\mu = \\mu_d + \\frac(\\mu_s-\\mu_d)(1+(\\frac(\\dot{\\epsilon}_{ii})(\\dot{\\epsilon}_C)))^x$  "
                           "where $\\mu_s$ and $\\mu_d$ are the friction angle at low and high strain rates, "
                           "respectively. $\\dot{\\epsilon}_{ii}$ is the second invariant of the strain rate and "
                           "$\\dot{\\epsilon}_C$ is the 'dynamic characteristic strain rate' where $\\mu = (\\mu_s+\\mu_d)/2$. "
                           "The 'dynamic friction smoothness exponent' x controls how "
                           "smooth or step-like the change from $\\mu_s$ to $\\mu_d$ is. "
                           "The equation is modified after Equation (13) in \\cite{van_dinther_seismic_2013}. "
                           "$\\mu_s$ and $\\mu_d$ can be specified by setting 'Angles of internal friction' and "
                           "'Dynamic angles of internal friction', respectively. "
                           "This relationship is similar to rate-and-state friction constitutive relationships, which "
                           "are applicable to the strength of rocks during earthquakes.");

        // Dynamic friction paramters
        prm.declare_entry ("Dynamic characteristic strain rate", "1e-12",
                           Patterns::Double (0),
                           "The characteristic strain rate value, where the angle of friction takes the middle "
                           "between the dynamic and the static angle of friction. When the effective strain rate "
                           "in a cell is very high, the dynamic angle of friction is taken, when it is very low "
                           "the static angle of internal friction is chosen. Around the dynamic characteristic "
                           "strain rate, there is a smooth gradient from the static to the dynamic friction "
                           "angle. "
                           "Units: \\si{\\per\\second}.");

        prm.declare_entry ("Dynamic angles of internal friction", "2",
                           Patterns::List(Patterns::Double(0)),
                           "List of dynamic angles of internal friction, $\\phi$, for background material and compositional "
                           "fields, for a total of N+1 values, where N is the number of compositional fields. "
                           "Dynamic angles of friction are used as the current friction angle when the effective "
                           "strain rate in a cell is well above the characteristic strain rate. "
                           "Units: \\si{\\degree}.");

        prm.declare_entry ("Dynamic friction smoothness exponent", "1",
                           Patterns::Double (0),
                           "An exponential factor in the equation for the calculation of the friction angle "
                           "when a static and a dynamic friction angle are specified. A factor of 1 returns the equation "
                           "to Equation (13) in \\cite{van_dinther_seismic_2013}. A factor between 0 and 1 makes the "
                           "curve of the friction angle vs. the strain rate more smooth, while a factor $>$ 1 makes "
                           "the change between static and dynamic friction angle more steplike. "
                           "Units: none.");
      }



      template <int dim>
      void
      FrictionOptions<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Get the number of fields for composition-dependent material properties
        // including the background field.
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        // Friction dependence parameters
        if (prm.get ("Friction dependence mechanism") == "none")
          friction_dependence_mechanism = independent;
        else if (prm.get ("Friction dependence mechanism") == "dynamic friction")
          friction_dependence_mechanism = dynamic_friction;
        else
          AssertThrow(false, ExcMessage("Not a valid friction dependence option!"));

        // Dynamic friction parameters
        dynamic_characteristic_strain_rate = prm.get_double("Dynamic characteristic strain rate");

        dynamic_angles_of_internal_friction = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Dynamic angles of internal friction"))),
                                                                                      n_fields,
                                                                                      "Dynamic angles of internal friction");
        // Convert angles from degrees to radians
        for (unsigned int i = 0; i<dynamic_angles_of_internal_friction.size(); ++i)
          {
            AssertThrow(dynamic_angles_of_internal_friction[i] <= 90, ExcMessage("Dynamic angles of friction must be <= 90 degrees"));
            dynamic_angles_of_internal_friction[i] *= numbers::PI/180.0;
          }

        dynamic_friction_smoothness_exponent = prm.get_double("Dynamic friction smoothness exponent");
      }
    }
  }
}



// Explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class FrictionOptions<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
