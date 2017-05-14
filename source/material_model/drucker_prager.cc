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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/material_model/drucker_prager.h>
#include <aspect/utilities.h>
#include <aspect/newton.h>

using namespace dealii;

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
      //set up additional output for the derivatives
      MaterialModel::MaterialModelDerivatives<dim> *derivatives;
      derivatives = out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim> >();

      /*if (derivatives == NULL)
      {

      const unsigned int n_points = out.viscosities.size();
      const unsigned int n_comp = out.reaction_terms[0].size();
      out.additional_outputs.push_back(
        std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
        (new MaterialModel::MaterialModelDerivatives<dim> (n_points)));
      }

      derivatives = out.template get_additional_output<MaterialModelDerivatives<dim> >();*/

      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // const Point<dim> position = in.position[i];
          const double temperature = in.temperature[i];
          const double pressure=std::max(in.pressure[i],0.0);

          // calculate effective viscosity
          if (in.strain_rate.size())
            {
              // This function calculates viscosities assuming that all the compositional fields
              // experience the same strain rate (isostrain). Since there is only one process in
              // this material model (a general powerlaw) we do not need to worry about how to
              // distribute the strain-rate and stress over the processes.
              //double composition_viscosities;
              SymmetricTensor<2,dim> composition_viscosity_strain_rate_derivatives;//(volume_fractions.size());
              double composition_viscosity_pressure_derivatives = 0;//(volume_fractions.size());

              const SymmetricTensor<2,dim> strain_rate_deviator = deviator(in.strain_rate[i]);
              const double edot_ii_strict = std::sqrt(0.5*strain_rate_deviator*strain_rate_deviator);

              const double sqrt3 = std::sqrt(3.0);
              const double degree_to_rad = numbers::PI/180;
              //const int ct = -1;

              //for (unsigned int c=0; c < volume_fractions.size(); ++c)
              //{
              // If strain rate is zero (like during the first time step) set it to some very small number
              // to prevent a division-by-zero, and a floating point exception.
              // Otherwise, calculate the square-root of the norm of the second invariant of the deviatoric-
              // strain rate (often simplified as epsilondot_ii)
              //const double edot_ii = 2.0 * std::max(edot_ii_strict, min_strain_rate[c] * min_strain_rate[c]);

              // Find effective viscosities for each of the individual phases
              // Viscosities should have same number of entries as compositional fields

              // Power law equation
              // edot_ii_i = A_i * stress_ii_i^{n_i} * d^{-m} \exp\left(-\frac{E_i^* + PV_i^*}{n_iRT}\right)
              // where ii indicates the square root of the second invariant and
              // i corresponds to diffusion or dislocation creep

              // The isostrain condition implies that the viscosity averaging should be arithmetic (see above).
              // We have given the user freedom to apply alternative bounds, because in diffusion-dominated
              // creep (where n_diff=1) viscosities are stress and strain-rate independent, so the calculation
              // of compositional field viscosities is consistent with any averaging scheme.
              // TODO: Mailed Bob and asked why the averaging should be arithmetic. Bob replyed that this has to
              //       do with effective medium theory. Have to look into this a bit more.
              //const double stress_exponent_inv = 1/stress_exponent[c];

              bool bool_min_strain_rate = false;
              double strain_rate_effective = edot_ii_strict;
              if (edot_ii_strict < reference_strain_rate * reference_strain_rate)
                {
                  strain_rate_effective = 2.0 * reference_strain_rate * reference_strain_rate;
                  bool_min_strain_rate = true;
                }

//std::cout << "i = " << i << "sr = " << strain_rate_effective << " -> " << in.strain_rate[0] << ", eta = " << eta_dislocation << ", A sei c = " << A_sei_c << ", exp_disloc = " << exp_dislocation << ", inner exp = " << (activation_energies_dislocation[c] + activation_volumes_dislocation[c] * in.pressure[i])*RT_inv*stress_exponent_inv << ", aed = " << activation_energies_dislocation[c] << ", avd = " << activation_volumes_dislocation[c] << ", p = " << in.pressure[i] << ", RT_inv = " << RT_inv <<  std::endl;
              // plasticity
              const double sin_phi = std::sin(angle_of_internal_friction * degree_to_rad);
              const double cos_phi = std::cos(angle_of_internal_friction * degree_to_rad);
              const double strain_rate_effective_inv = 1/strain_rate_effective;
              const double strength_inv_part = 6.0/(sqrt3*(3.0-sin_phi));
              const double strength = (dim == 2 ?
                                       pressure * sin_phi + cohesion * cos_phi
                                       :
                                       (cohesion*cos_phi*strength_inv_part)
                                       + (sin_phi * strength_inv_part) * pressure);
              const double eta_plastic     = 0.5 * strength * strain_rate_effective_inv;
              //std::cout << "eta plastic = " << eta_plastic << " = 0.5 * " << strength << " * " << strain_rate_effective_inv << std::endl;
              //std::cout << "strenght  = " << strength << " = (" << cohesions[c] << "*" << cos_phi << "*" << strength_inv_part << ") + " << (sin_phi * strength_inv_part) * pressure << std::endl;


              /*
              if (dim==3)
              {
              const int ct = -1;
              strength = ((6.0*cohesion*std::cos(phi))/(std::sqrt(3.0)*(3.0+ct*std::sin(phi)))) +
                 ((6.0*std::sin(phi))/(std::sqrt(3.0)*(3.0+ct*std::sin(phi)))) * std::max(pressure,0.0);
              }
              else strength = std::max(pressure,0.0) * std::sin(phi) + cohesion * std::cos(phi);
              return strength / (2.0*strain_rate_norm);
               */


              const double composition_viscosity = std::max(std::min(eta_plastic, maximum_viscosity), minimum_viscosity);

              Assert(dealii::numbers::is_finite(composition_viscosity),ExcMessage ("Error: Viscosity is not finite."));
              if (/*this->get_parameters().newton_theta != 0 &&*/ derivatives != NULL)
                {
                  //std::cout << "viscosity = " << composition_viscosities[c] << ", min_visc = " << min_visc[c] << ", viscoplast = " << std::min(eta_viscous,eta_plastic) << ", visco = " << eta_viscous << ", plastic = " << eta_plastic << ", eta_diffusion = " << eta_diffusion << ", eta_dislocation = " << eta_dislocation << std::endl;
                  //std::cout << "viscosity = " << composition_viscosities[c] << ", min_visc = " << min_visc[c]  << ", max_visc = " << max_visc[c] << ", bool = " << bool_min_strain_rate << std::endl;
                  if (!bool_min_strain_rate && composition_viscosity < maximum_viscosity && composition_viscosity > minimum_viscosity)
                    {
                      //std::cout << "viscous = " << eta_viscous << ", eta_plastic = " << eta_plastic << std::endl;

                      //std::cout << "I am fully plastic! " << strain_rate_deviator  << ", " << in.strain_rate[i] << std::endl;
                      composition_viscosity_strain_rate_derivatives = 0.5 * (eta_plastic * (-1/(edot_ii_strict * edot_ii_strict))) * strain_rate_deviator;
                      composition_viscosity_pressure_derivatives = (dim == 2 ?
                                                                    0.5 * sin_phi * strain_rate_effective_inv
                                                                    :
                                                                    (sin_phi*strength_inv_part));
                      /*std::cout << "B " << eta_viscous << ":" << eta_plastic << ", eIIstict = " << edot_ii_strict << "; "
                      << -1/(edot_ii_strict * edot_ii_strict) << ", " << 0.5 * (eta_plastic * (-1/(edot_ii_strict * edot_ii_strict))) *
                       strain_rate_deviator  << " ; " << sin_phi << "*" << strength_inv_part << ", sr = "
                       << composition_viscosities_strain_rate_derivatives[c] << std::endl;*/
                    }
                  else
                    {
                      composition_viscosity_strain_rate_derivatives = 0;
                      //std::cout << "Set " << i << " to 0, becasue of bool_min_strain_rate = " << bool_min_strain_rate << " = " << edot_ii_strict << "<" << min_strain_rate[c] * min_strain_rate[c] << ", max:" << max_visc[c] << ", min:" << min_visc[c]  << std::endl;
                    }
                  //std::cout << std::setprecision (15) << "B " << composition_viscosities[c] << " = " << eta_viscous << ":" << eta_plastic << ", " << composition_viscosities_strain_rate_derivatives[c] << " = " << sin_phi << "*" << strength_inv_part << std::endl;
                }

              //}

              out.viscosities[i] = composition_viscosity;//Utilities::weighted_p_norm_average(volume_fractions, composition_viscosities, viscosity_averaging_p);
              /*std::cout << "C " << out.viscosities[i] << " = ";
              for(unsigned int c=0; c < volume_fractions.size(); ++c)
                std::cout << ", " << composition_viscosities[c];
              std::cout << std::endl;*/
              Assert(dealii::numbers::is_finite(out.viscosities[i]),ExcMessage ("Error: Averaged viscosity is not finite."));

              if (/*this->get_parameters().newton_theta != 0 &&*/ derivatives != NULL)
                {
                  //std::cout << "A derivatives->dviscosities_dpressure[i] = " << derivatives->dviscosities_dpressure[i] ;
                  /*for(unsigned int c=0; c < volume_fractions.size(); ++c)
                    std::cout << ", " << composition_viscosities_pressure_derivatives[c];
                  std::cout << std::endl;*/
                  derivatives->viscosity_derivative_wrt_strain_rate[i] = composition_viscosity_strain_rate_derivatives;//Utilities::derivatives_weighed_p_norm_average(out.viscosities[i],volume_fractions, composition_viscosities, composition_viscosities_strain_rate_derivatives, viscosity_averaging_p);

                  derivatives->viscosity_derivative_wrt_pressure[i] = composition_viscosity_pressure_derivatives;//Utilities::derivatives_weighed_p_norm_average(out.viscosities[i],volume_fractions, composition_viscosities, composition_viscosities_pressure_derivatives, viscosity_averaging_p);
//p_norm_average(volume_fractions, composition_viscosities_pressure_derivatives, viscosity_averaging_p);
                  //std::cout << "B derivatives->dviscosities_dpressure[i] = " << derivatives->dviscosities_dpressure[i] << ", " << composition_viscosities_pressure_derivatives[0] << ":" << composition_viscosities_pressure_derivatives[1] << ":" << composition_viscosities_pressure_derivatives[2] << ":" << composition_viscosities_pressure_derivatives[3] << ":" << composition_viscosities_pressure_derivatives[4] << std::endl;
                  //std::cout << derivatives->dviscosities_dpressure[i] << ":" << composition_viscosities_pressure_derivatives[0] << ":" << std::endl;
                  /*#ifdef DEBUG
                                    for (int x = 0; x < dim; x++)
                                      for (int y = 0; y < dim; y++)
                                        if (!dealii::numbers::is_finite(derivatives->dviscosities_dstrain_rate[i][x][y]))
                                          std::cout << "Error: Averaged viscosity to strain-rate devrivative is not finite." << std::endl;

                                    //derivatives->dviscosities_dpressure[i]    = 0;* /
                                    if (!dealii::numbers::is_finite(derivatives->dviscosities_dpressure[i]))
                                      {
                                        std::cout << "Error: Averaged viscosity to pressure devrivative is not finite. " << std::endl;
                                        for (unsigned int c=0; c < volume_fractions.size(); ++c)
                                          std::cout << composition_viscosities_pressure_derivatives[c] << ",";
                                        std::cout << std::endl;
                                      }
                                    Assert(dealii::numbers::is_finite(derivatives->dviscosities_dpressure[i]),ExcMessage ("Error: Averaged dviscosities_dpressure is not finite."));
                                    for (int x = 0; x < dim; x++)
                                      for (int y = 0; y < dim; y++)
                                        Assert(dealii::numbers::is_finite(derivatives->dviscosities_dstrain_rate[i][x][y]),ExcMessage ("Error: Averaged dviscosities_dstrain_rate is not finite."));
                                    //Assert(dealii::numbers::is_nan(out.viscosities[i]),ExcMessage ("Error: Averaged viscosity is not finite."));
                  #endif*/
                  /*if (isnan(derivatives->dviscosities_dpressure[i]))
                  {
                                      std::cout << "Error: Averaged viscosity to pressure devrivative is nam: " << derivatives->dviscosities_dpressure[i] << ":" << composition_viscosities_pressure_derivatives[0]  << ":" << composition_viscosities_pressure_derivatives[1] << ":" << composition_viscosities_pressure_derivatives[2] << ":" << composition_viscosities_pressure_derivatives[3] << ":" << composition_viscosities_pressure_derivatives[4]  << ":" << composition_viscosities_pressure_derivatives[5] << ":" << composition_viscosities_pressure_derivatives[6]  << ":" << composition_viscosities_pressure_derivatives[7] << ":" << composition_viscosities_pressure_derivatives[8]  << ":" << composition_viscosities_pressure_derivatives[9] << std::endl;
                                    AssertThrow(false,ExcMessage ("Error: Averaged viscosity is nan."));
                  }
                  if (!dealii::numbers::is_finite(derivatives->dviscosities_dpressure[i]))
                  {
                                      std::cout << "Error: Averaged viscosity to pressure devrivative is not finite: " << derivatives->dviscosities_dpressure[i] << ":" << composition_viscosities_pressure_derivatives[0] << std::endl;
                                    AssertThrow(false,ExcMessage ("Error: Averaged viscosity is not finite."));
                  }*/
                  /*#ifdef DEBUG
                                    for (int x = 0; x < dim; x++)
                                      for (int y = 0; y < dim; y++)
                                        if (!dealii::numbers::is_finite(derivatives->dviscosities_dstrain_rate[i][x][y]))
                                          std::cout << "Error: Averaged viscosity to strain-rate devrivative is not finite." << std::endl;

                                    //derivatives->dviscosities_dpressure[i]    = 0;
                                    / *if (isnan(derivatives->dviscosities_dpressure[i]))
                                      std::cout << "Error: Averaged viscosity to pressure devrivative is not finite: " << derivatives->dviscosities_dpressure[i] << ":" << composition_viscosities_pressure_derivatives[0] << std::endl;
                                    AssertThrow(false,ExcMessage ("Error: Averaged viscosity is not finite."))* /
                  #endif*/
                }
            }
          out.densities[i] = reference_rho * (1 - thermal_expansivity * (temperature - reference_T));
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
    viscosity (const double /*temperature*/,
               const double pressure,
               const std::vector<double> &/*composition*/,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &/*position*/) const
    {
      // For the very first time this function is called
      // (the first iteration of the first timestep), this function is called
      // with a zero input strain rate. We provide a representative reference
      // strain rate for this case, which avoids division by zero and produces
      // a representative first guess of the viscosities.
      // In later iterations and timesteps we calculate the second moment
      // invariant of the deviatoric strain rate tensor.
      // This is equal to the negative of the second principle
      // invariant calculated with the function second_invariant.
      const double strain_rate_dev_inv2 = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
                                            ?
                                            reference_strain_rate * reference_strain_rate
                                            :
                                            std::fabs(second_invariant(deviator(strain_rate))));

      // In later timesteps, we still need to care about cases of very small
      // strain rates. We expect the viscosity to approach the maximum_viscosity
      // in these cases. This check prevents a division-by-zero.
      if (std::sqrt(strain_rate_dev_inv2) <= std::numeric_limits<double>::min())
        return maximum_viscosity;

      // To avoid negative yield strengths and eventually viscosities,
      // we make sure the pressure is not negative
      const double strength = ( (dim==3)
                                ?
                                ( 6.0 * cohesion * std::cos(angle_of_internal_friction) + 2.0 * std::max(pressure,0.0) * std::sin(angle_of_internal_friction) )
                                / ( std::sqrt(3.0) * ( 3.0 + std::sin(angle_of_internal_friction) ) )
                                :
                                cohesion * std::cos(angle_of_internal_friction) + std::max(pressure,0.0) * std::sin(angle_of_internal_friction) );

      // Rescale the viscosity back onto the yield surface
      const double viscosity = strength / ( 2.0 * std::sqrt(strain_rate_dev_inv2) );

      // Cut off the viscosity between a minimum and maximum value to avoid
      // a numerically unfavourable large viscosity range.
      const double effective_viscosity = 1.0 / ( ( 1.0 / ( viscosity + minimum_viscosity ) ) + ( 1.0 / maximum_viscosity ) );

      return effective_viscosity;

    }


    template <int dim>
    double
    DruckerPrager<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }


    template <int dim>
    double
    DruckerPrager<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return reference_specific_heat;
    }


    template <int dim>
    double
    DruckerPrager<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return thermal_conductivities;
    }


    template <int dim>
    double
    DruckerPrager<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &) const
    {
      return reference_rho * (1 - thermal_expansivity * (temperature - reference_T));
    }


    template <int dim>
    double
    DruckerPrager<dim>::
    thermal_expansion_coefficient (const double,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return thermal_expansivity;
    }


    template <int dim>
    double
    DruckerPrager<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
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
          thermal_expansivity              = prm.get_double ("Thermal expansion coefficient");
          prm.enter_subsection ("Viscosity");
          {
            minimum_viscosity        = prm.get_double ("Minimum viscosity");
            maximum_viscosity        = prm.get_double ("Maximum viscosity");
            reference_strain_rate    = prm.get_double ("Reference strain rate");
            // Convert degrees to radians
            angle_of_internal_friction = prm.get_double ("Angle of internal friction") * numbers::PI/180.0;
            cohesion                 = prm.get_double ("Cohesion");
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
                                   "Section~\\ref{parameters:Material_20model/Drucker_20Prager}."
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
