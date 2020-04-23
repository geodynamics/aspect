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

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/newton.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * This is the material model used for the viscoplastic convection
     * benchmark models that were presented in the open access Tosi et al. (2015) paper:
     * @code
     *  @Article{T15,
     *   Author = {Tosi, N. and Stein, C. and Noack, L. and H\"uttig, C. and Maierova, P. and Samual, H. and Davies, D. R. and Wilson, C. R. and Kramer, S. C. and Thieulot, C. and Glerum, A. and Fraters, M. and Spakman, W. and Rozel, A. and Tackley, P. J.},
     *    Title = {A community benchmark for viscoplastic thermal convection in a 2-D square box},
     *    Journal = {Geochemistry, Geophysics, Geosystems},
     *    Year = {2015}
     *    Pages = {2175-2196},
     *    Volume = {16},
     *    Doi = {10.1002/2015GC005807}}
     * @endcode
     *
     * It is a material model that consists of globally constant values for all
     * material parameters except density and viscosity. Density is temperature dependent,
     * but practically constant through a small value of the thermal expansivity,
     * and viscosity depends on temperature, depth and strain rate.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class NewtonTosiMaterial : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        virtual
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * Method to declare parameters related to depth-dependent model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to depth-dependent model
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * Method that indicates whether material is compressible. Depth dependent model is compressible
         * if and only if base model is compressible.
         */
        virtual bool is_compressible () const;

        /**
         * Method to calculate reference viscosity for the depth-dependent model. The reference
         * viscosity is determined by evaluating the depth-dependent part of the viscosity at
         * the mean depth of the model.
         */
        virtual double reference_viscosity () const;



      private:

        double viscosity (const double                  temperature,
                          const double                  pressure,
                          const std::vector<double>    &compositional_fields,
                          const SymmetricTensor<2,dim> &strain_rate,
                          const Point<dim>             &position) const;

        double density (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position) const;
        /*
         * Function to compute the linear viscosity
         * according to equation (7) of Tosi et al. 2015.
         */
        double viscolin(const double etaT,
                        const double etaZ,
                        const double T,
                        const double depth) const;

        /*
         * Function to compute the plastic viscosity
         * according to equation (8) of the paper.
         */
        double viscoplast(const double eta_astrix,
                          const double stress_y,
                          const double strain_rate_norm) const;

        bool use_analytical_derivative;

        /**
         * The density at reference temperature
         */
        double reference_rho;

        /*
         * The reference temperature
         */
        double reference_T;

        /*
         * The thermal expansivity
         */
        double thermal_alpha;

        /*
         * The reference specific heat
         */
        double reference_specific_heat;

        /**
         * The thermal conductivity.
         */
        double thermal_k;


        /**
         * As described in Tosi et al (2015), the viscosity \eta is computed as the
         * harmonic average of a linear and nonlinear part.
         *
         * The linear part is calculated as follows
         * (see below for the meaning of the used parameter names):
         * $\eta_{lin}(T,z) = \exp(-\ln(\text{eta\_T} * T + \ln(\text{eta\_Z}) * z)$
         * while the strain rate dependent nonlinear part is computed as:
         * $\eta_{plast}(\dot\epsilon) = \text{eta\_astrix} + \frac{\text{sigma\_yield}}{\sqrt(\dot\epsilon:\dot\epsilon)}$
         */

        /*
         * The reference viscosity
         */
        double eta;

        /*
         * The linear viscosity parameter pertaining to
         * the viscosity contrast due to temperature
         */
        double eta_T;

        /*
         * The linear viscosity parameter pertaining to
         * the viscosity contrast due to pressure (depth)
         */
        double eta_Z;

        /*
         * The effective viscosity at high stresses that is
         * part of the plastic viscosity
         */
        double eta_astrix;

        /*
         * The constant yield stress that is used in the
         * plastic viscosity formula
         */
        double sigma_yield;

        /*
         * The lower viscosity cut-off value
         */
        double eta_minimum;

        /*
         * The upper viscosity cut-off value
         */
        double eta_maximum;

        /*
         * The initial guess of the viscosity
         */
        double eta_initial;


    };

  }
}



#include <aspect/utilities.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    NewtonTosiMaterial<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      //set up additional output for the derivatives
      MaterialModelDerivatives<dim> *derivatives;
      derivatives = out.template get_additional_output<MaterialModelDerivatives<dim> >();
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          //const double depth = this->get_geometry_model().depth(in.position[i]);
          // as documented, if the strain rate array is empty, then do not compute the
          // viscosities
          if (in.strain_rate.size() > 0)
            out.viscosities[i]                  = viscosity                     (in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);

          out.densities[i]                      = density                       (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.thermal_expansion_coefficients[i] = thermal_alpha;
          out.specific_heat[i]                  = reference_specific_heat;
          out.thermal_conductivities[i]         = thermal_k;
          out.compressibilities[i]              = 0;

          if (derivatives != NULL && in.strain_rate.size() > 0)
            {
              if (use_analytical_derivative)
                {
                  if ((/*this->get_timestep_number() == 0 &&*/ in.strain_rate[i].norm() == 0))
                    {
                      derivatives->viscosity_derivative_wrt_strain_rate[i] = 0;
                      derivatives->viscosity_derivative_wrt_pressure[i] = 0;
                    }
                  else
                    {
                      const double deviator_strain_rate_norm = in.strain_rate[i].norm();//deviator(in.strain_rate[i]).norm();
                      const double part1 = (eta_astrix * deviator_strain_rate_norm + sigma_yield);

                      derivatives->viscosity_derivative_wrt_strain_rate[i] = -(0.5*out.viscosities[i]*out.viscosities[i]) * (sigma_yield/(part1 * part1 * deviator_strain_rate_norm)) * in.strain_rate[i];
                      derivatives->viscosity_derivative_wrt_pressure[i] = 0;
                    }
                }
              else
                {
                  // finite difference
                  const double finite_difference_accuracy = 1e-7;
                  SymmetricTensor<2,dim> zerozero = SymmetricTensor<2,dim>();
                  SymmetricTensor<2,dim> onezero = SymmetricTensor<2,dim>();
                  SymmetricTensor<2,dim> oneone = SymmetricTensor<2,dim>();

                  zerozero[0][0] = 1;
                  onezero[1][0]  = 0.5; // because symmetry doubles this entry
                  oneone[1][1]   = 1;

                  SymmetricTensor<2,dim> strain_rate_zero_zero = in.strain_rate[i] + std::fabs(in.strain_rate[i][0][0]) * finite_difference_accuracy * zerozero;
                  SymmetricTensor<2,dim> strain_rate_one_zero = in.strain_rate[i] + std::fabs(in.strain_rate[i][1][0]) * finite_difference_accuracy * onezero;
                  SymmetricTensor<2,dim> strain_rate_one_one = in.strain_rate[i] + std::fabs(in.strain_rate[i][1][1]) * finite_difference_accuracy * oneone;

                  double eta_zero_zero = viscosity(in.temperature[i], in.pressure[i], in.composition[i], strain_rate_zero_zero, in.position[i]);
                  double deriv_zero_zero = eta_zero_zero - out.viscosities[i];

                  if (deriv_zero_zero != 0)
                    {
                      if (strain_rate_zero_zero[0][0] != 0)
                        {
                          deriv_zero_zero /= std::fabs(strain_rate_zero_zero[0][0]) * finite_difference_accuracy;
                        }
                      else
                        {
                          deriv_zero_zero = 0;
                        }

                    }

                  double eta_one_zero = viscosity(in.temperature[i], in.pressure[i], in.composition[i], strain_rate_one_zero, in.position[i]);
                  double deriv_one_zero = eta_one_zero - out.viscosities[i];

                  if (deriv_one_zero != 0)
                    {
                      if (strain_rate_one_zero[1][0] != 0)
                        {
                          deriv_one_zero /= std::fabs(strain_rate_one_zero[1][0]) * finite_difference_accuracy;
                        }
                      else
                        {
                          deriv_one_zero = 0;
                        }
                    }

                  double eta_one_one = viscosity(in.temperature[i], in.pressure[i], in.composition[i], strain_rate_one_one, in.position[i]);
                  double deriv_one_one = eta_one_one - out.viscosities[i];

                  if (eta_one_one != 0)
                    {
                      if (strain_rate_one_one[1][1] != 0)
                        {
                          deriv_one_one /= std::fabs(strain_rate_one_one[1][1]) * finite_difference_accuracy;
                        }
                      else
                        {
                          deriv_one_one = 0;
                        }
                    }

                  derivatives->viscosity_derivative_wrt_strain_rate[i][0][0] = deriv_zero_zero;
                  derivatives->viscosity_derivative_wrt_strain_rate[i][1][0] = deriv_one_zero;
                  derivatives->viscosity_derivative_wrt_strain_rate[i][1][1] = deriv_one_one;


                  derivatives->viscosity_derivative_wrt_pressure[i] = 0;//deriv_pressure;

                }
            }
        }
    }
    /*
     * Function to compute the linear viscosity
     * according to equation (7) of Tosi et al. 2015.
     */
    template <int dim>
    double
    NewtonTosiMaterial<dim>::
    viscolin(const double etaT,
             const double etaZ,
             const double T,
             const double z) const
    {
      return std::exp((-1.0 * std::log(etaT) * T ) + (std::log(etaZ) * z));
    }

    /*
     * Function to compute the plastic viscosity
     * according to equation (8) of the paper.
     */
    template <int dim>
    double
    NewtonTosiMaterial<dim>::
    viscoplast(const double etaastrix,
               const double stressy,
               const double strainratenorm) const
    {
      return etaastrix +(stressy/strainratenorm);
    }

    /*
     * Function to calculate the viscosity according
     * to equation (6) of the paper.
     */
    template <int dim>
    double
    NewtonTosiMaterial<dim>::
    viscosity (const double temperature,
               const double,
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &) const
    {

      // In the first nonlinear iteration of the (pre-refinement steps of the) first time step,
      // strain rate is zero, so we set viscosity to eta_initial, a user-defined guess of the viscosity.
      if (strain_rate.norm() == 0)
        {
          return eta_initial;
        }

      // Otherwise we compute the linear viscosity and, if needed, the plastic viscosity.
      double viscosity = 0.0;
      Point<dim> test();
      const double visc_linear = viscolin(eta_T,eta_Z,temperature,0);

      if (eta_astrix == 0.0 && sigma_yield == 0.0)
        {
          viscosity = visc_linear;
        }
      else
        {
          const double visc_plastic = viscoplast(eta_astrix,sigma_yield,strain_rate.norm());//deviator(strain_rate).norm());

          // Compute the harmonic average (equation (6) of the paper)
          viscosity = 2.0 / ((1.0 / visc_linear) + (1.0 / visc_plastic));
        }

      // Cut-off the viscosity by user-defined values to avoid possible very large viscosity ratios
      viscosity = std::max(std::min(viscosity,eta_maximum),eta_minimum);

      return viscosity;
    }


    template <int dim>
    double
    NewtonTosiMaterial<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    NewtonTosiMaterial<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &,
             const Point<dim> &) const
    {
      return reference_rho * (1.0 - thermal_alpha * (temperature - reference_T));
    }



    template <int dim>
    bool
    NewtonTosiMaterial<dim>::
    is_compressible () const
    {
      return false;
    }


    template <int dim>
    void
    NewtonTosiMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      // Default values are for Case 1 of Tosi et al. (2015).
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Newton Tosi model");
        {
          prm.declare_entry ("Use analytical derivative", "false",
                             Patterns::Bool (),
                             "Wether to use the analytical or the finite difference derivative for the Newton method.");
          prm.declare_entry ("Reference density", "1",
                             Patterns::Double (0),
                             "The value of the reference density $\\rho_0$.");
          prm.declare_entry ("Reference temperature", "0",
                             Patterns::Double (0),
                             "The value of the reference temperature $T_0$. The reference temperature is used "
                             "in the density calculation.");
          prm.declare_entry ("Reference viscosity", "1e-1",
                             Patterns::Double (0),
                             "The value of the constant reference viscosity $\\eta_0$.");
          prm.declare_entry ("Minimum viscosity", "1e-6",
                             Patterns::Double (0),
                             "The value of the minimum cut-off viscosity $\\eta_min$.");
          prm.declare_entry ("Maximum viscosity", "1e1",
                             Patterns::Double (0),
                             "The value of the maximum cut-off viscosity $\\eta_max$.");
          prm.declare_entry ("Initial viscosity", "1e-1",
                             Patterns::Double (0),
                             "The value of the initial viscosity guess $\\eta_init$.");
          prm.declare_entry ("Thermal conductivity", "1",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$.");
          prm.declare_entry ("Reference specific heat", "1",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$.");
          prm.declare_entry ("Thermal expansion coefficient", "1e-6",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\alpha$.");
          prm.declare_entry ("Thermal viscosity parameter", "1e5",
                             Patterns::Double (0),
                             "The value of the thermal viscosity parameter $\\eta_T$, "
                             "as used in equation (7) of the paper.");
          prm.declare_entry ("Pressure viscosity parameter", "1e0",
                             Patterns::Double (0),
                             "The value of the pressure viscosity parameter $\\eta_Z$, "
                             "as used in equation (7) of the paper.");
          prm.declare_entry ("Yield stress", "0",
                             Patterns::Double (0),
                             "The value of the plastic viscosity yield stress $\\sigma_yield$, "
                             "as used in equation (8) of the paper.");
          prm.declare_entry ("Nonlinear viscosity constant", "0",
                             Patterns::Double (0),
                             "The value of the plastic viscosity constant $\\eta_astrix$, "
                             "as used in equation (8) of the paper.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    NewtonTosiMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Newton Tosi model");
        {
          use_analytical_derivative  = prm.get_bool ("Use analytical derivative");
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          eta                        = prm.get_double ("Reference viscosity");
          thermal_k                  = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          eta_T                      = prm.get_double ("Thermal viscosity parameter");
          eta_Z                      = prm.get_double ("Pressure viscosity parameter");
          sigma_yield                = prm.get_double ("Yield stress");
          eta_astrix                 = prm.get_double ("Nonlinear viscosity constant");
          eta_minimum                = prm.get_double ("Minimum viscosity");
          eta_maximum                = prm.get_double ("Maximum viscosity");
          eta_initial                = prm.get_double ("Initial viscosity");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(NewtonTosiMaterial,
                                   "Newton Tosi material",
                                   "This is the same as the Tosi becmark material, but has the derivatives needed for the "
                                   "Newton solver. "
                                   "A material model that has constant values "
                                   "for all coefficients but the density and viscosity as described in the "
                                   "open access paper of Tosi et al. 2015. "
                                   "The default parameter values are chosen according to Case 1 of the paper. "
                                   "All of the values that define this model are read "
                                   "from a section ``Material model/Newton Tosi model'' in the input file, see "
                                   "Section~\\ref{parameters:Material_model/Newton_Tosi_model}."
                                   "\n\n"
                                   "This model uses the following set of equations for the two coefficients that "
                                   "are non-constant (see equation (6) - (10) of the paper): "
                                   "\\begin{align}"
                                   "  \\eta(T,z,\\dot \\epsilon) &= 2\\frac{1}{\\frac{1}{\\eta_{lin}(T,z)}\\frac{1}{\\eta_{plast}(\\dot\\epsilon}}, \\\\"
                                   "  \\rho(T) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0,"
                                   "\\end{align}"
                                   "where z represents depth."
                                   "\n\n"
                                   "The linear and plastic viscosity parts are defined as follows:"
                                   "\\begin{align}"
                                   "  \\eta_{lin}(T,z) &= \\exp(-\\ln(\\eta_T)T+\\ln(\\eta_z)z), \\\\"
                                   "  \\eta_{plast}(\\dot\\epsilon) &= \\eta^{*}+\\frac{\\sigma_y}{\\sqrt(\\dot\\epsilon : \\dot\\epsilon)} "
                                   "\\end{align} "
                                   "\n\n"
                                   "Note that this model uses the formulation that assumes an incompressible "
                                   "medium despite the fact that the density follows the law "
                                   "$\\rho(T)=\\rho_0(1-\\alpha(T-T_0))$. ")
  }
}
