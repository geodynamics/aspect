/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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

#include <aspect/global.h>
#include <aspect/initial_temperature/anelasticity_vs_to_temperature.h>
#include <aspect/material_model/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>

#include <boost/lexical_cast.hpp>
#include <boost/math/tools/minima.hpp>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    AnelasticVsToTemperature<dim>::AnelasticVsToTemperature ()
    {}

    template <int dim>
    void
    AnelasticVsToTemperature<dim>::initialize ()
    {
      Utilities::AsciiDataInitial<dim>::initialize(1);
    }

    template <int dim>
    double
    AnelasticVsToTemperature<dim>::
    fVs (const double x,
         const double      depth,
         const double      absolute_Vs,
         const double      mu0,
         const double      dmudT,
         const double      dmudP,
         const double      viscosity_prefactor,
         const double      activation_energy,
         const double      activation_volume,
         const double      solidus_gradient,
         const bool      use_original_model) const
    {
      return std::abs(yamauchi_takei_Vs(x,depth,mu0,dmudT,dmudP,viscosity_prefactor,activation_energy,activation_volume
                                        ,solidus_gradient,use_original_model)-absolute_Vs);
    }

    template <int dim>
    double
    AnelasticVsToTemperature<dim>::
    fdV (const double x, const double bulk_modulus, const double bulk_modulus_pressure_derivative, const double pressure ) const
    {
      return std::abs((bulk_modulus*(3./2.)*(std::pow(x,7./3.)-std::pow(x,5./3.))*(1+(((3./4.)*
                                                                                      (bulk_modulus_pressure_derivative-4))*(std::pow(x,2./3.)-1))))-pressure);
    }

    // set up initial temperature
    template <int dim>
    double
    AnelasticVsToTemperature<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // determine depth
      const double depth = this->get_geometry_model().depth(position);

      // declare temperature
      double temperature;

      // read absolute Vs in from ascii file
      const double absolute_Vs = Utilities::AsciiDataInitial<dim>::get_data_component(position,0);

      if (depth >= no_perturbation_depth)
        {
          // convert absolute Vs into temperature
          // check if using Yamauchi & Takei 2016 parameterization
          // specify Brent algorithm parameters
          const double a=273;
          const double b=3273;
          // create fVs function to use in Brent minimization and calculate temperature
          auto bfunc = [ &,this] (const double x)
          {
            return fVs(x, depth, absolute_Vs, mu0, dmudT,dmudP,
                       viscosity_prefactor,activation_energy,activation_volume,solidus_gradient,use_original_model);
          };
          // determine maximum Vs
          const double maximum_Vs=yamauchi_takei_Vs(273.,depth,mu0,dmudT,dmudP,viscosity_prefactor,activation_energy,activation_volume
                                                    ,solidus_gradient,use_original_model);

          // where absolute Vs exceeds maximum Vs, set temperature to 273 K
          if (absolute_Vs>maximum_Vs)
            {
              temperature=273.;
            }
          else
            {
              temperature=boost::math::tools::brent_find_minima(bfunc,a,b,16).first;
            }
        }
      else
        {
          // set temperature to constant above specified depth
          temperature = upper_temperature;
        }
      // return the absolute temperature in Kelvin
      return temperature;
    }

    // set up function to calculate shear wave velocity from Yamauchi & Takei 2016 parameterisation
    template <int dim>
    double
    AnelasticVsToTemperature<dim>::
    yamauchi_takei_Vs (const double temperature,
                       const double depth,
                       const double mu0,
                       const double dmudT,
                       const double dmudP,
                       const double viscosity_prefactor,
                       const double activation_energy,
                       const double activation_volume,
                       const double solidus_gradient,
                       const bool use_original_model) const
    {
      // specify anelasticity parameters
      const double critical_homologous_temperature = 0.94;
      const double reduction_factor = 5;
      const double background_amplitude = 0.664;
      const double background_slope = 0.38;
      const double peak_period = 6e-5;
      const double melt_viscosity_factor = 0;
      const double melt_peak_factor = 0;
      const double reference_temperature = 1473;
      const double reference_pressure = 1.5e9;
      const double grain_size = 1e-3;
      const double reference_grain_size = 1e-3;
      const double grain_size_exponent = 3;
      const double pressure_gradient = 3e-5;
      const double gas_constant=8.3145;
      // specify Grose & Afonso (2013) density parameters
      const double a=1;
      const double b=3;
      const double bulk_modulus=130e9;
      const double bulk_modulus_pressure_derivative=4.8;
      const double gruneisen_parameter=6;
      const double reference_density=3330;
      // specify original density parameters
      const double original_density=3291;
      const double original_thermal_expansivity=3.59e-5;
      const double original_bulk_modulus=115.2;
      // initialize solidus
      const double T_solidus = 1326.0 + 273 + (((depth-50000)/1e3)*solidus_gradient);
      // initialize homologous_temperature
      const double homologous_temperature = temperature/T_solidus;
      // initialize pressures
      const double pressure = depth/pressure_gradient;
      // declare other parameters
      double viscosity,viscosity_reduction_factor,unrelaxed_compliance;
      // begin calculation of Vs
      if (homologous_temperature<critical_homologous_temperature)
        {
          viscosity_reduction_factor=1.0;
        }
      else if ((homologous_temperature>=critical_homologous_temperature) && (homologous_temperature<1))
        {
          viscosity_reduction_factor=std::exp((-1.0*((homologous_temperature-critical_homologous_temperature)/(homologous_temperature-
                                                     (homologous_temperature*critical_homologous_temperature))))*std::log(reduction_factor));
        }
      else
        {
          viscosity_reduction_factor=(1.0/reduction_factor)*std::exp(-melt_viscosity_factor);
        }
      viscosity = std::pow(grain_size/reference_grain_size,grain_size_exponent)*viscosity_prefactor*std::exp((activation_energy/gas_constant)
                  *(1.0/temperature-1.0/reference_temperature))*std::exp((activation_volume/gas_constant)*(pressure/temperature-reference_pressure/
                                                                         reference_temperature))*viscosity_reduction_factor;
      unrelaxed_compliance=1.0/(1e9*(mu0+(dmudP*pressure*1e-9)+(dmudT*(temperature-273))));

      if (temperature<273)
        {
          // Vs is too high to give realistic temperature so viscosity, attenuation and unrelaxed compliance are reset
          viscosity=1e40;
          unrelaxed_compliance=1./(1e9*(mu0+(dmudP*pressure*1e-9)));
        }
      // evaluate Maxwell normalised shear wave period
      const double maxwell_relaxation_time=viscosity*unrelaxed_compliance;
      double period;
      if (use_original_model == true)
        {
          // set shear wave period as constant
          period=100;
        }
      else
        {
          // calculate shear wave period incorporating depth dependence of Forsyth 1992
          period=(3*depth)/4200;
        }
      // determine normalised shear wave period
      const double normalised_period=period/(2*M_PI*maxwell_relaxation_time);
      // determine peak amplitudes
      double peak_amplitude;
      if (homologous_temperature < 0.91)
        {
          peak_amplitude=0.01;
        }
      else if ((homologous_temperature >= 0.91) && (homologous_temperature < 0.96))
        {
          peak_amplitude=0.01+(0.4*(homologous_temperature-0.91));
        }
      else if ((homologous_temperature >= 0.96) && (homologous_temperature < 1))
        {
          peak_amplitude=0.03;
        }
      else
        {
          peak_amplitude=0.03+melt_peak_factor;
        }
      // determine peak widths
      double peak_width;
      if (homologous_temperature < 0.92)
        {
          peak_width=4;
        }
      else if ((homologous_temperature >= 0.92) && (homologous_temperature < 1))
        {
          peak_width=4+(37.5*(homologous_temperature-0.92));
        }
      else
        {
          peak_width=7;
        }
      // determine density
      double compressibility,pressure_dependent_density,integrated_thermal_expansivity,density,isothermal_volume_change;
      if (use_original_model == true)
        {
          // calculate density using original parameters
          density=original_density*(1-(original_thermal_expansivity*((temperature-273)-600))+((pressure*1e-9)/original_bulk_modulus));
        }
      else
        {
          // create fdV function to use in Brent minimization and calculate isothermal_volume_change and density using
          // expressions in Grose & Afonso 2013
          auto vfunc = [ &,this] (const double x)
          {
            return fdV(x, bulk_modulus, bulk_modulus_pressure_derivative, pressure);
          };
          isothermal_volume_change=boost::math::tools::brent_find_minima(vfunc,a,b,16).first;
          compressibility=isothermal_volume_change*std::exp((gruneisen_parameter+1)*(std::pow(isothermal_volume_change,-1)-1));
          pressure_dependent_density=reference_density*isothermal_volume_change;
          integrated_thermal_expansivity=(2.832e-5*(temperature-273))+((0.758e-8/2)*(std::pow(temperature,2)-std::pow(273,2)));
          density=pressure_dependent_density*(1-(compressibility*integrated_thermal_expansivity));
        }
      double storage_compliance,anelastic_Vs;
      // determine J1 term (real part of complex compliance)
      storage_compliance=unrelaxed_compliance*(1+((background_amplitude*std::pow(normalised_period,background_slope))
                                                  /background_slope)+((std::sqrt(2*M_PI)/2)*peak_amplitude*peak_width*(1-
                                                                      std::erf((std::log(peak_period/normalised_period))/(std::sqrt(2)*peak_width)))));
      // calculate Vs
      anelastic_Vs=1/(std::sqrt(density*storage_compliance)*1e3);
      return anelastic_Vs;
    }

    template <int dim>
    void
    AnelasticVsToTemperature<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.declare_entry ("Remove temperature heterogeneity down to specified depth",
                           boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                           Patterns::Double (),
                           "This will remove density and viscosity variations down to the specified depth "
                           "(in meters). Note that your resolution has to be adequate to capture this cutoff. "
                           "For example if you specify a depth of 660km, but your closest spherical depth layers "
                           "are only at 500km and 750km (due to a coarse resolution) it will only zero out "
                           "heterogeneities down to 500km. Similar caution has to be taken when using adaptive meshing.");
        prm.declare_entry ("Use original density and frequency model of Yamauchi and Takei", "true",
                           Patterns::Bool(),
                           "Use original density and frequency model of Yamauchi and Takei where density"
                           "has simple pressure-dependence and shear wave is period set at 100s");
        prm.declare_entry ("Upper temperature", "600",
                           Patterns::Double(),
                           "The value of temperature in the above the zero perturbation depth. Units: $K$.");
        prm.declare_entry ("Reference shear modulus", "72.45",
                           Patterns::Double(),
                           "The value of the unrelaxed shear modulus at surface pressure-temperature conditions $\\mu_U^0$. Units: $GPa$.");
        prm.declare_entry ("Temperature derivative of shear modulus", "-0.01094",
                           Patterns::Double(),
                           "The value of the temperature derivative of shear modulus $\\dmu_U/dT$. Units: $GPa K^{-1}$.");
        prm.declare_entry ("Pressure derivative of shear modulus", "1.987",
                           Patterns::Double(),
                           "The value of the pressure derivative of shear modulus $\\dmu_U/dP$.");
        prm.declare_entry ("Viscosity prefactor", "6.22e21",
                           Patterns::Double(),
                           "The value of the viscosity prefactor $\\eta_r$. Units: $Pa s$.");
        prm.declare_entry ("Activation energy", "462.5e3",
                           Patterns::Double(),
                           "The value of the activation energy $\\E_a$. Units: $J/mol$.");
        prm.declare_entry ("Activation volume", "7.913e-6",
                           Patterns::Double(),
                           "The value of the activation volume $\\V_a$. Units: $m^3/mol$.");
        prm.declare_entry ("Solidus gradient", "1.018",
                           Patterns::Double(),
                           "The value of the solidus gradient $\\ $dT_s/dz}$. Units: $K/km$.");


        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
                                                          "box_2d_Vs_YT16.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    AnelasticVsToTemperature<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        no_perturbation_depth   = prm.get_double ("Remove temperature heterogeneity down to specified depth");
        use_original_model = prm.get_bool ("Use original density and frequency model of Yamauchi and Takei");
        upper_temperature = prm.get_double("Upper temperature");
        mu0 = prm.get_double ("Reference shear modulus");
        dmudT = prm.get_double ("Temperature derivative of shear modulus");
        dmudP = prm.get_double ("Pressure derivative of shear modulus");
        viscosity_prefactor = prm.get_double ("Viscosity prefactor");
        activation_energy = prm.get_double ("Activation energy");
        activation_volume = prm.get_double ("Activation volume");
        solidus_gradient = prm.get_double ("Solidus gradient");

        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(AnelasticVsToTemperature,
                                              "anelasticity Vs to temperature",
                                              "Implementation of a model in which the initial temperature is "
                                              "calculated from files containing absolute shear wave velocity (Vs) "
                                              "data in ascii format. This plug-in allows you to convert Vs into anelastic"
                                              "temperature, accounting for the behaviour of mantle material. "
                                              "Currently only the Yamauchi & Takei 2016 parameterization is implemented. "
                                              "Note the required format of the "
                                              "input data: The first lines may contain any number of comments "
                                              "if they begin with `#', but one of these lines needs to "
                                              "contain the number of grid points in each dimension as "
                                              "for example `# POINTS: 3 3'. "
                                              "The order of the data columns "
                                              "has to be `x', `y', `Velocity  [km/s]' in a 2d model and "
                                              " `x', `y', `z', `Velocity [km/s]' in a 3d model, which means that "
                                              "there has to be a single column "
                                              "containing the temperature. "
                                              "Note that the data in the input "
                                              "files need to be sorted in a specific order: "
                                              "the first coordinate needs to ascend first, "
                                              "followed by the second and the third at last in order to "
                                              "assign the correct data to the prescribed coordinates. "
                                              "If you use a spherical model, "
                                              "then the assumed grid changes. `x' will be replaced by "
                                              "the radial distance of the point to the bottom of the model, "
                                              "`y' by the azimuth angle and `z' by the polar angle measured "
                                              "positive from the north pole. The grid will be assumed to be "
                                              "a latitude-longitude grid. Note that the order "
                                              "of spherical coordinates is `r', `phi', `theta' "
                                              "and not `r', `theta', `phi', since this allows "
                                              "for dimension independent expressions.")
  }
}
