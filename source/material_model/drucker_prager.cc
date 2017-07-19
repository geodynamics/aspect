/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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


#include <aspect/material_model/drucker_prager.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>


namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    DruckerPrager<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // set up additional output for the derivatives
      MaterialModel::MaterialModelDerivatives<dim> *derivatives;
      derivatives = out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const double temperature = in.temperature[i];
          // To avoid negative yield strengths and eventually viscosities,
          // we make sure the pressure is not negative
          const double pressure=std::max(in.pressure[i],0.0);

          // calculate effective viscosity
          if (in.strain_rate.size() > 0)
            {
              const SymmetricTensor<2,dim> strain_rate_deviator = deviator(in.strain_rate[i]);

              // For the very first time this function is called
              // (the first iteration of the first timestep), this function is called
              // with a zero input strain rate. We provide a representative reference
              // strain rate for this case, which avoids division by zero and produces
              // a representative first guess of the viscosities.
              //
              // In later iterations and timesteps we calculate the second moment
              // invariant of the deviatoric strain rate tensor.
              // This is equal to the negative of the second principle
              // invariant of the deviatoric strain rate (calculated with the function second_invariant),
              // as shown in Appendix A of Zienkiewicz and Taylor (Solid Mechanics, 2000).
              //
              // The negative of the second principle invariant is equal to 0.5 e_dot_dev_ij e_dot_dev_ji,
              // where e_dot_dev is the deviatoric strain rate tensor. The square root of this quantity
              // gives the common definition of effective strain rate.
              const double edot_ii_strict = (this->simulator_is_initialized() == false
                                             ?
                                             // no simulator object available -- we are probably in a unit test
                                             std::fabs(second_invariant(strain_rate_deviator))
                                             :
                                             // simulator object is available, but we need to treat the
                                             // first time step separately
                                             ((this->get_timestep_number() == 0)
                                              &&
                                              (in.strain_rate[i].norm() <= std::numeric_limits<double>::min())
                                              ?
                                              reference_strain_rate * reference_strain_rate
                                              :
                                              std::fabs(second_invariant(strain_rate_deviator))));

              const double sqrt3 = std::sqrt(3.0);

              const double strain_rate_effective = edot_ii_strict;

              // In later time steps, we still need to care about cases of very small
              // strain rates. We expect the viscosity to approach the maximum_viscosity
              // in these cases. This check prevents a division-by-zero.
              if (std::sqrt(strain_rate_effective) <= std::numeric_limits<double>::min())
                {
                  out.viscosities[i] = maximum_viscosity;

                  if (derivatives != NULL)
                    {
                      derivatives->viscosity_derivative_wrt_strain_rate[i] = 0.0;
                      derivatives->viscosity_derivative_wrt_pressure[i] = 0.0;
                    }
                }
              else
                {
                  // plasticity
                  const double sin_phi = std::sin(angle_of_internal_friction);
                  const double cos_phi = std::cos(angle_of_internal_friction);
                  const double strain_rate_effective_inv = 1./(2.*std::sqrt(strain_rate_effective));
                  const double strength_inv_part = 1./(sqrt3*(3.0+sin_phi));

                  const double strength = ( (dim==3)
                                            ?
                                            ( 6.0 * cohesion * cos_phi + 2.0 * pressure * sin_phi) * strength_inv_part
                                            :
                                            cohesion * cos_phi + pressure * sin_phi );

                  // Rescale the viscosity back onto the yield surface
                  const double eta_plastic = strength * strain_rate_effective_inv;

                  // Cut off the viscosity between a minimum and maximum value to avoid
                  // a numerically unfavourable large viscosity range.
                  const double effective_viscosity = 1.0 / ( ( 1.0 / ( eta_plastic + minimum_viscosity ) ) + ( 1.0 / maximum_viscosity ) );

                  Assert(dealii::numbers::is_finite(effective_viscosity), ExcMessage ("Error: Viscosity is not finite."));

                  out.viscosities[i] = effective_viscosity;

                  Assert(dealii::numbers::is_finite(out.viscosities[i]),
                         ExcMessage ("Error: Averaged viscosity is not finite."));

                  if (derivatives != NULL)
                    {
                      const double averaging_factor = effective_viscosity * effective_viscosity * std::pow(eta_plastic + minimum_viscosity,-2);
                      const SymmetricTensor<2,dim> effective_viscosity_strain_rate_derivatives
                        = averaging_factor * 0.5 * (eta_plastic * (-1.0/std::sqrt(edot_ii_strict * edot_ii_strict))) * strain_rate_deviator;
                      const double effective_viscosity_pressure_derivatives = averaging_factor * (dim == 2 ?
                                                                                                  sin_phi * strain_rate_effective_inv
                                                                                                  :
                                                                                                  (sin_phi*strength_inv_part));

                      derivatives->viscosity_derivative_wrt_strain_rate[i] = effective_viscosity_strain_rate_derivatives;

                      derivatives->viscosity_derivative_wrt_pressure[i] = effective_viscosity_pressure_derivatives;

                      Assert(dealii::numbers::is_finite(derivatives->viscosity_derivative_wrt_pressure[i]),
                             ExcMessage ("Error: Averaged viscosity_derivative_wrt_pressure is not finite."));
                      for (int x = 0; x < dim; x++)
                        for (int y = 0; y < dim; y++)
                          Assert(dealii::numbers::is_finite(derivatives->viscosity_derivative_wrt_strain_rate[i][x][y]),
                                 ExcMessage ("Error: Averaged viscosity_derivative_wrt_strain_rate is not finite."));

                    }
                }
            }
          out.densities[i] = reference_rho * (1.0 - thermal_expansivity * (temperature - reference_T));
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          // Specific heat at the given positions.
          out.specific_heat[i] = reference_specific_heat;
          // Thermal conductivity at the given positions.
          out.thermal_conductivities[i] = thermal_conductivities;
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
          for (unsigned int c=0; c < in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }



    template <int dim>
    double
    DruckerPrager<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    bool
    DruckerPrager<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    DruckerPrager<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Drucker Prager");
        {
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "The reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in the density calculation. Units: $K$.");
          prm.declare_entry ("Reference viscosity", "1e22",
                             Patterns::Double (0),
                             "The value of the reference viscosity $\\eta_0$. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Reference specific heat", "1250",
                             Patterns::Double (0),
                             "The value of the specific heat $C_p$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
          prm.enter_subsection ("Viscosity");
          {
            prm.declare_entry ("Minimum viscosity", "1e19",
                               Patterns::Double (0),
                               "The value of the minimum viscosity cutoff $\\eta_min$. Units: $Pa\\;s$.");
            prm.declare_entry ("Maximum viscosity", "1e24",
                               Patterns::Double (0),
                               "The value of the maximum viscosity cutoff $\\eta_max$. Units: $Pa\\;s$.");
            prm.declare_entry ("Reference strain rate", "1e-15",
                               Patterns::Double (0),
                               "The value of the initial strain rate prescribed during the "
                               "first nonlinear iteration $\\dot{\\epsilon}_ref$. Units: $1/s$.");
            prm.declare_entry ("Angle of internal friction", "0",
                               Patterns::Double (0),
                               "The value of the angle of internal friction $\\phi$. "
                               "For a value of zero, in 2D the von Mises "
                               "criterion is retrieved. Angles higher than 30 degrees are "
                               "harder to solve numerically. Units: degrees.");
            prm.declare_entry ("Cohesion", "2e7",
                               Patterns::Double (0),
                               "The value of the cohesion $C$. Units: $Pa$.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    DruckerPrager<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Drucker Prager");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          reference_eta              = prm.get_double ("Reference viscosity");
          thermal_conductivities     = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_expansivity        = prm.get_double ("Thermal expansion coefficient");
          prm.enter_subsection ("Viscosity");
          {
            minimum_viscosity          = prm.get_double ("Minimum viscosity");
            maximum_viscosity          = prm.get_double ("Maximum viscosity");
            reference_strain_rate      = prm.get_double ("Reference strain rate");
            // Convert degrees to radians
            angle_of_internal_friction = prm.get_double ("Angle of internal friction") * numbers::PI/180.0;
            cohesion                   = prm.get_double ("Cohesion");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.viscosity = NonlinearDependence::strain_rate;
      this->model_dependence.density = NonlinearDependence::none;

      if (angle_of_internal_friction==0.0)
        this->model_dependence.viscosity |= NonlinearDependence::pressure;

      if (thermal_expansivity != 0)
        this->model_dependence.density = NonlinearDependence::temperature;
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(DruckerPrager,
                                   "drucker prager",
                                   "A material model that has constant values "
                                   "for all coefficients but the density and viscosity. The defaults for all "
                                   "coefficients are chosen to be similar to what is believed to be correct "
                                   "for Earth's mantle. All of the values that define this model are read "
                                   "from a section ``Material model/Drucker Prager'' in the input file, see "
                                   "Section~\\ref{parameters:Material_20model/Drucker_20Prager}. "
                                   "Note that the model does not take into account any dependencies of "
                                   "material properties on compositional fields. "
                                   "\n\n"
                                   "The viscosity is computed according to the Drucker Prager frictional "
                                   "plasticity criterion (non-associative) based on a user-defined "
                                   "internal friction angle $\\phi$ and cohesion $C$. In 3D: "
                                   " $\\sigma_y = \\frac{6 C \\cos(\\phi)}{\\sqrt(3) (3+\\sin(\\phi))} + "
                                   "\\frac{2 P \\sin(\\phi)}{\\sqrt(3) (3+\\sin(\\phi))}$, "
                                   "where $P$ is the pressure. "
                                   "See for example Zienkiewicz, O. C., Humpheson, C. and Lewis, R. W. (1975), "
                                   "G\\'{e}otechnique 25, No. 4, 671-689. "
                                   "With this formulation we circumscribe instead of inscribe the Mohr Coulomb "
                                   "yield surface. "
                                   "In 2D the Drucker Prager yield surface is the same "
                                   "as the Mohr Coulomb surface: "
                                   " $\\sigma_y = P \\sin(\\phi) + C \\cos(\\phi)$. "
                                   "Note that in 2D for $\\phi=0$, these criteria "
                                   "revert to the von Mises criterion (no pressure dependence). "
                                   "See for example Thieulot, C. (2011), PEPI 188, 47-68. "
                                   "\n\n"
                                   "Note that we enforce the pressure to be positive to prevent negative "
                                   "yield strengths and viscosities. "
                                   "\n\n"
                                   "We then use the computed yield strength to scale back the viscosity on "
                                   "to the yield surface using the Viscosity Rescaling Method described in "
                                   "Kachanov, L. M. (2004), Fundamentals of the Theory of Plasticity, "
                                   "Dover Publications, Inc. (Not Radial Return.)"
                                   "A similar implementation can be found in GALE "
                                   "(https://geodynamics.org/cig/software/gale/gale-manual.pdf). "
                                   "\n\n"
                                   "To avoid numerically unfavourably large (or even negative) viscosity ranges, "
                                   "we cut off the viscosity with a user-defined minimum and maximum viscosity: "
                                   "$\\eta_eff = \\frac{1}{\\frac{1}{\\eta_min + \\eta}+ "
                                   "\\frac{1}{\\eta_max}}$. "
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\beta(T-T_{\\text{ref}}))$. ")
  }
}
