/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: billen_2013.cc Sept-Oct, 2013 Katrina M Arredondo $  */


#include <aspect/material_model/billen_2013.h>
#include <deal.II/base/parameter_handler.h>
#include <iostream>
#include <algorithm>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    double
    Billen_2013<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &composition,       /*composition*/
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &position) const
    {
		const double delta_temp = temperature-reference_T;
		double viscosity;
		const double Re =  6371.0e3; //meters
		const double R =  8.3145;  //%J/mol K
       		const double depth = Re - position.norm();
		const double y = Re * atan(position(0)/position(1));
 
        const double F_diff = 1e0/(std::pow(2e0, ((stress_exponent_diffusion-1e0)/stress_exponent_diffusion))*(std::pow(3e0, ((stress_exponent_diffusion+1e0)/(2*stress_exponent_diffusion)))));
        const double F_dis = 1e0/(std::pow(2e0, ((stress_exponent_dislocation-1e0)/stress_exponent_dislocation))*(std::pow(3e0, ((stress_exponent_dislocation+1e0)/(2*stress_exponent_dislocation)))));

        const double betaA = 4.3e-12; // 1/Pa adiabatic compressibility (Burkett and Billen , 2009)
        const double Plc_660 = (-1e0/betaA)*std::log(1e0 - reference_rho*g_0*betaA*660e3); // Pa - Turcotte book eq. 4-250
        const double Plc_0 = (-1e0/betaA)*std::log(1e0 - reference_rho*g_0*betaA*0e0); // Pa - Turcotte book eq. 4-250
        const double Pm5 = (Plc_660 - Plc_0)/(660e3 - 0e0); // identical to CitcomS configuration
        double eedot, tau, tau_new, eta_p, eta_new, m_ad, viscosity_diffusion;
      
//	double nonnewt = std::sqrt(velocity_values.norm() * velocity_values.norm(); 
	// CALCULATE EFFECTIVE STRAIN RATE 
	double tmst = this->get_timestep_number();	

        //if(tmst > 0){	
        	/* get second invariant for all elements */
		if (dim == 2)
            		eedot = std::pow(0.5 * (std::pow(strain_rate[0][0], 2) + std::pow(strain_rate[1][1], 2)) + std::pow(strain_rate[0][1], 2), 0.5);
		else{
			std::cout << "This model currently only supports 2D models." << endl;
			eedot = 1.0e-15;
		}
        //}
        //else{
            /* initialize with unity if no velocities around */
        //    eedot = 1.0e-15; // 1/s 		
            //fprintf(stderr,"Starting with Fixed Strain-rate or 1e-15 s^1: sdepv_visit = %d\n", tmst);
            //std::cout << "Starting with Fixed Strain-rate or 1e-15 s^1: timestep = " << tmst << endl;
        //}
        
        if (adiabat_add == true)
            m_ad = 3e-4; // = 0.3 K/km = 0.3/1000 K/m; Fixed mantle temperature gradient 
        else
            m_ad = 0e0;
		
        // ********************************************************** //
		// Compute viscosity
		
		// default value
        viscosity = eta_0;
		
		// Ask for temperature, stress and pressure dependence.
        if (TDEPV == true)
        {
            if (depth <= 660.0e3)
            {
                const double visc_diff_one = F_diff * std::pow(std::pow(grain_size_diffusion_UM , grain_size_exponent_diffusion) / (prefactor_diffusion * std::pow(water_content_diffusion , water_exponent_diffusion)) , (1e0 / stress_exponent_diffusion));
                const double visc_diff_two = (activation_energy_diffusion + depth*Pm5*activation_volume_diffusion_UM)/(stress_exponent_diffusion*R*(temperature + m_ad * depth));
                viscosity_diffusion = visc_diff_one * std::exp(visc_diff_two);
	//	if ((depth > 100e3) && (depth < 120e3))
	//		std::cout << "Hello, from the viscosity file." << " eedot: " << eedot << endl;
		//	std::cout << "depth: " << depth << " temperature: " << temperature << " m_ad: " << m_ad << " m_ad*depth: " << m_ad*depth << "visc_diff_one: " << visc_diff_one << " visc_diff_two: " << visc_diff_two << " viscosity_diffusion: " << viscosity_diffusion << endl;
            }
            if (depth > 660.0e3)
            {
                const double visc_diff_one_LM = F_diff * std::pow(std::pow(grain_size_diffusion_LM , grain_size_exponent_diffusion) / (prefactor_diffusion * std::pow(water_content_diffusion , water_exponent_diffusion)) , (1e0 / stress_exponent_diffusion));
                const double visc_diff_two_LM = (activation_energy_diffusion + depth*Pm5*activation_volume_diffusion_LM)/(stress_exponent_diffusion*R*(temperature + m_ad * depth));
                viscosity_diffusion = visc_diff_one_LM * std::exp(visc_diff_two_LM);
		//if ((depth > 660e3) && (depth < 710e3))
		//	std::cout << "depth: " << depth << " temperature: " << temperature << " visc_diff_one_LM: " << visc_diff_one_LM << " visc_diff_two_LM: " << visc_diff_two_LM << " viscosity_diffusion: " << viscosity_diffusion << endl;
            }
            
            viscosity = std::min(1e50, viscosity_diffusion);
            //if ((depth > 100e3) && (depth < 150e3))
		//std::cout << "Viscosity Newtonian: " << viscosity << " at the depth: " << depth << " with viscosity_diffusion: " << viscosity_diffusion << " temperature: " << temperature << "\n";
	}
            
        if ((SDEPV == true) && (depth <= 660.0e3))
        {
            const double visc_dis_one = F_dis * std::pow(1e0/(prefactor_dislocation * std::pow(water_content_dislocation, water_exponent_dislocation)) , (1e0/stress_exponent_dislocation));
            const double visc_dis_two = (activation_energy_dislocation + depth*Pm5*activation_volume_dislocation_UM)/(stress_exponent_dislocation*R*(temperature + m_ad * depth));
            const double viscosity_dislocation = visc_dis_one * std::pow(eedot,((1e0-stress_exponent_dislocation)/stress_exponent_dislocation)) * std::exp(visc_dis_two);
                
            const double etads = std::min(1e50, viscosity_dislocation);
            const double etadf = viscosity;
            
            viscosity = etadf*etads/(etadf+etads);
        }
            
        if (PDEPV == true)
        {
            // depth dependent yield stress
            tau = min_yield_stress + depth*rate_yield_stress;
            tau_new = std::min(tau, max_yield_stress);
                
            // viscosity
            eta_p = tau_new/(2.0*eedot);
            eta_new = std::min(viscosity, eta_p);
	    //if ((depth > 10e3) && (depth < 30e3) && (y > 4000e3) && (y < 4050e3))
            //    std::cout << "Prior Viscosity: " << viscosity << " at the depth: " << depth << " tau: " << tau << " tau_new: " << tau_new << " max_yield_stress : " << max_yield_stress << " min_yield_stress: " << min_yield_stress << " eta_p " << eta_p << " eedot: " << eedot << " eta_new: " << eta_new << "\n"; 
            viscosity = eta_new;
        }
        
        // If there is a weakzone and kinematic boundary condition at the surface, need to add condition for top nodes
        // ADD THIS LATERi
        //if ((depth > 100e3) && (depth < 150e3))
        //        std::cout << "Viscosity Newtonian: " << viscosity << " at the depth: " << depth << " temperature: " << temperature << "\n";
	
	viscosity = std::min(max_visc, viscosity);
	const double right_extent = trench_location + weakzone_depth/tan(weakzone_angle * numbers::PI / 180) + weakzone_width/2;

	if (weakzone == true)
	{
		if ((y < right_extent) && (y > (trench_location - (weakzone_width/2))))
		{
			if (depth < weakzone_depth)
			{
				const double weakzone_center = trench_location + depth/tan(weakzone_angle * numbers::PI / 180); // calculate weakzone center for a defined angle from the surface
				const double y_weakzone = y - weakzone_center; 		// distance from the center of the weakzone
				const double ts = 0.5*(weakzone_width - weakest_width); 	// distance over which viscosity decreases
				const double tp = 0.5 * (weakest_width + ts); 	   
				const double y_temp = std::abs(y_weakzone);
				const double dn = (y_temp - tp + 0.5*ts) / ts; 		// nondimensional distance
				double S;

				// Define center of weakzone and the sides
				if (std::abs(y_weakzone) < (tp - 0.5*ts))
					S = 1;
				else{
					S = 1 - std::pow(dn,2)*(3 - 2*dn);	// function for weakzone viscosity
					if (std::abs(y_weakzone) > weakzone_width/2)
						S = 0;
				}

				// S should be between 0 and 1 
				//std::cout << "S: " << S << " dn: " << dn << " ts: " << ts << " y_weakzone: " << y_weakzone << " weakzone_center: " << weakzone_center << " y: " << y << " depth: " << depth << " weakzone_angle: " << weakzone_angle << std::endl;
				//const double etatemp = eta_0 * std::pow(10, (std::log10(viscosity / eta_0) * (1 - S)));
				//std::cout << "Weakzone Viscosity" << etatemp << ". Nondim Visc: " << std::log10(viscosity / eta_0) << std::endl;
				const double etatemp = std::pow(10, (std::log10(viscosity) - S * (std::log10(viscosity) - std::log10(weakzone_viscosity))));
				viscosity = std::min(viscosity, etatemp);
			}
		}
	}
		
	viscosity = std::max(min_visc, viscosity);
	
	return viscosity;
    }


    template <int dim>
    double
    Billen_2013<dim>::
    reference_viscosity () const
    {
      return eta_0;
    }

    template <int dim>
    double
    Billen_2013<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    Billen_2013<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return alpha_0;
    }

    template <int dim>
    double
    Billen_2013<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return cp_0;
    }

    template <int dim>
    double
    Billen_2013<dim>::
    reference_cp () const
    {
      return cp_0;
    }

    template <int dim>
    double
    Billen_2013<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return K_0*(reference_rho*cp_0);
    }

    template <int dim>
    double
    Billen_2013<dim>::
    reference_thermal_diffusivity () const
    {
        return K_0;
    }

    template <int dim>
    double
    Billen_2013<dim>::
    density (const double temperature,
             const double,
             const std::vector<double> &compositional_fields, /*composition*/
             const Point<dim> &) const
    {
      return (reference_rho * (1 - alpha_0 * (temperature - reference_T)));
    }


    template <int dim>
    double
    Billen_2013<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return alpha_0;
    }


    template <int dim>
    double
    Billen_2013<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }

    template <int dim>
    bool
    Billen_2013<dim>::
    viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the viscosity() function
      // to see the dependencies
        return ((dependence & NonlinearDependence::pressure)
                ||
                (dependence & NonlinearDependence::temperature)
                ||
                (dependence & NonlinearDependence::strain_rate));
    }

    template <int dim>
    bool
    Billen_2013<dim>::
    density_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      // compare this with the implementation of the density() function
      // to see the dependencies
      if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
          &&
          (alpha_0 != 0))
        return true;
      else
        return false;
    }

    template <int dim>
    bool
    Billen_2013<dim>::
    compressibility_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Billen_2013<dim>::
    specific_heat_depends_on (const NonlinearDependence::Dependence) const
    {
      return false;
    }

    template <int dim>
    bool
    Billen_2013<dim>::
    thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
    {
      return false;
    }


    template <int dim>
    bool
    Billen_2013<dim>::
    is_compressible () const
    {
      return false;
    }



    template <int dim>
    void
    Billen_2013<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Billen_2013 model");
        {
			prm.declare_entry ("Temperature Dependent", "true",
							   Patterns::Bool (),
							   "Temperature Dependence");
			prm.declare_entry ("Stress Dependent", "true",
							   Patterns::Bool (),
							   "Stress Dependence");
			prm.declare_entry ("Pressure Dependent", "true",
							   Patterns::Bool (),
							   "Pressure Dependence");
			prm.declare_entry ("Constant mantle temperature", "true",
							   Patterns::Bool (),
							   "Constant mantle temperature? true");
                        prm.declare_entry ("Include Weakzone", "false",
                                                           Patterns::Bool (),
                                                           "Include a weakzone");
          prm.declare_entry ("Weakzone Viscosity", "10^21",
                             Patterns::Double (),
                             "Viscosity of the weakzone. Units: $Pa s$.");
          prm.declare_entry ("Weakzone Angle", "30e3",
                             Patterns::Double (0),
                             "Angle of weakzone from the surface, defined in radians. Units: $radians$.");
          prm.declare_entry ("Weakzone Width", "30e3",
                             Patterns::Double (0),
                             "Weakzone width defined in meters. Units: $m$.");
         prm.declare_entry ("Trench Location", "3000e3",
                             Patterns::Double (0),
                             "Location of weakzone center at the surface defined in meters from the left grid edge. Units: $m$.");
          prm.declare_entry ("Weakest Weakzone Width", "5e3",
                             Patterns::Double (0),
                             "Width of the weakest center of the weakzone defined in meters. Units: $m$.");
          prm.declare_entry ("Maximum Weakzone Depth", "80e3",
                             Patterns::Double (0),
                             "Maximum extent of the weakzone with depth defined in meters. Units: $m$.");

          prm.declare_entry ("Maximum Viscosity", "1e24",
                             Patterns::Double (0),
                             "Maximum Viscosity. Units: $Pa s$.");
          prm.declare_entry ("Minimum Viscosity", "1e18",
                             Patterns::Double (0),
                             "Minimum Viscosity. Units: $Pa s$.");

          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "1400",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. Units: $K$.");
          prm.declare_entry ("Reference Viscosity", "1e20",
                             Patterns::Double (0),
                             "The value of the constant viscosity. Units: $kg/m/s$.");
          prm.declare_entry ("Reference Thermal diffusivity", "1e-6",
                             Patterns::Double (0),
                             "The value of the thermal diffusivity $K$. "
                             "Units: $m^2/s$.");
          prm.declare_entry ("Reference Heat Capacity", "1250",
                             Patterns::Double (0),
                             "The value of the heat capacity $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("Reference thermal expansion coefficient", "2e-5",
                             Patterns::Double (0),
                             "The value of the thermal expansion coefficient $\\beta$. "
                             "Units: $1/K$.");
		  prm.declare_entry ("Reference gravity", "9.81",
							   Patterns::Double (0),
							   "The value of the gravity constant."
							   "Units: $m/s^2$.");
			
			prm.declare_entry ("Activation energy diffusion", "317",
							   Patterns::Double (0),
							   "Activation energy diffusion $E$. Units: $kJ/mol$.");
			prm.declare_entry ("Activation volume diffusion in upper mantle", "4.0e-6",
							   Patterns::Double (0),
							   "Activation volume diffusion in upper mantle $V$. Units: $m^3/mol$.");
			prm.declare_entry ("Activation volume diffusion in lower mantle", "1.5e-6",
							   Patterns::Double (0),
							   "Activation volume diffusion in lower mantle $V$. Units: $m^3/mol$.");
			prm.declare_entry ("Stress exponent diffusion", "1e0",
							   Patterns::Double (0),
							   "Stress exponent diffusion $n$. Units: $none$.");
			prm.declare_entry ("Prefactor diffusion", "1e0",
							   Patterns::Double (0),
							   "Prefactor diffusion $A$. Units: $1/sPa^n\etam^pC_OH^-r$.");
			prm.declare_entry ("Grain Size Exponent diffusion", "3e0",
							   Patterns::Double (0),
							   "Grain Size Exponent diffusion $p$. Units: $none$.");
			prm.declare_entry ("Grain Size diffusion in upper mantle", "10000",
							   Patterns::Double (0),
							   "Grain Size diffusion in upper mantle $d$. Units: $micrometer$.");
			prm.declare_entry ("Grain Size diffusion in lower mantle", "110.4e3",
							   Patterns::Double (0),
							   "Grain Size diffusion in lower mantle $d$. Units: $micrometer$.");
			prm.declare_entry ("Water content diffusion", "1000",
							   Patterns::Double (0),
							   "Water content diffusion $C_OH$. Units: $ppm H/Si$.");
			prm.declare_entry ("Water exponent diffusion", "1e0",
							   Patterns::Double (0),
							   "Water exponent diffusion $r$. Units: $none$.");
			
			prm.declare_entry ("Activation energy dislocation", "496",
							   Patterns::Double (0),
							   "Activation energy dislocation $E$. Units: $kJ/mol$.");
			prm.declare_entry ("Activation volume dislocation in upper mantle", "11e-6",
							   Patterns::Double (0),
							   "Activation volume dislocation in upper mantle $V$. Units: $m^3/mol$.");
			prm.declare_entry ("Stress exponent dislocation", "3.5",
							   Patterns::Double (0),
							   "Stress exponent dislocation $n$. Units: $none$.");
			prm.declare_entry ("Prefactor dislocation", "9.0e-20",
							   Patterns::Double (0),
							   "Prefactor dislocation $A$. Units: $1/(sPa^n\etam^pC_OH^-r)$.");
			prm.declare_entry ("Water content dislocation", "1000",
							   Patterns::Double (0),
							   "Water content dislocation $C_OH$. Units: $ppm H/Si$.");
			prm.declare_entry ("Water exponent dislocation", "1.2",
							   Patterns::Double (0),
							   "Water exponent dislocation $r$. Units: $none$.");
			
            prm.declare_entry ("Initial Yield Stress", "1000e6",
							   Patterns::Double (0),
							   "Initial Yield Stress $\\sigma_a$. Units: $Pa$.");
            prm.declare_entry ("Yield Stress Change Rate", "1000e6",
							   Patterns::Double (0),
							   "Yield Stress Change Rate $\\sigma_b$. Units: $Pa/m$.");
			prm.declare_entry ("Maximum Yield Stress", "1000e6",
							   Patterns::Double (0),
							   "Maximum Yield Stress $\\sigma_y$. Units: $Pa$.");
		}
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Billen_2013<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Billen_2013 model");
        {
		TDEPV                 = prm.get_bool ("Temperature Dependent");
		SDEPV                 = prm.get_bool ("Stress Dependent");
		PDEPV                 = prm.get_bool ("Pressure Dependent");
		adiabat_add 	      = prm.get_bool ("Constant mantle temperature");
		weakzone              = prm.get_bool ("Include Weakzone");
		weakzone_viscosity        = prm.get_double ("Weakzone Viscosity");
		weakzone_angle        = prm.get_double ("Weakzone Angle");
		weakzone_width        = prm.get_double ("Weakzone Width");
		trench_location	      = prm.get_double ("Trench Location");
                weakest_width         = prm.get_double ("Weakest Weakzone Width");
		weakzone_depth        = prm.get_double ("Maximum Weakzone Depth");

		max_visc              = prm.get_double ("Maximum Viscosity");
		min_visc              = prm.get_double ("Minimum Viscosity");	
	
         	reference_rho              = prm.get_double ("Reference density");
         	reference_T                = prm.get_double ("Reference temperature");
          	eta_0                       = prm.get_double ("Reference Viscosity");
          	K_0						= prm.get_double ("Reference Thermal diffusivity");
          	cp_0						= prm.get_double ("Reference Heat Capacity");
          	alpha_0					= prm.get_double ("Reference thermal expansion coefficient");
          	g_0						= prm.get_double ("Reference gravity");;
			
		activation_energy_diffusion   = prm.get_double ("Activation energy diffusion");
		activation_volume_diffusion_UM   = prm.get_double ("Activation volume diffusion in upper mantle");
		activation_volume_diffusion_LM   = prm.get_double ("Activation volume diffusion in lower mantle");
		stress_exponent_diffusion          = prm.get_double ("Stress exponent diffusion");
		prefactor_diffusion					= prm.get_double ("Prefactor diffusion");
		grain_size_diffusion_UM				= prm.get_double ("Grain Size diffusion in upper mantle");
		grain_size_diffusion_LM				= prm.get_double ("Grain Size diffusion in lower mantle");
		grain_size_exponent_diffusion		= prm.get_double ("Grain Size Exponent diffusion");
		water_content_diffusion				= prm.get_double ("Water content diffusion");
		water_exponent_diffusion			= prm.get_double ("Water exponent diffusion");
			
			activation_energy_dislocation		= prm.get_double ("Activation energy dislocation");
			activation_volume_dislocation_UM	= prm.get_double ("Activation volume dislocation in upper mantle");
			stress_exponent_dislocation			= prm.get_double ("Stress exponent dislocation");
			prefactor_dislocation				= prm.get_double ("Prefactor dislocation");
			water_content_dislocation			= prm.get_double ("Water content dislocation");
			water_exponent_dislocation			= prm.get_double ("Water exponent dislocation");
			
            min_yield_stress = prm.get_double ("Initial Yield Stress");
            rate_yield_stress = prm.get_double ("Yield Stress Change Rate");
			max_yield_stress = prm.get_double ("Maximum Yield Stress");

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
    ASPECT_REGISTER_MATERIAL_MODEL(Billen_2013,
                                   "billen_2013",
                                   "A simple material model that has constant values "
                                   "for all coefficients but the density and viscosity. "
                                   "This model uses an incompressible medium."
                                   "The viscosity is dependent on temperature, pressure and"
                                   " stress, including a yield stress."
                                   "The value for the components of this formula and additional "
                                   "parameters are read from the paper 'Billen and Hirth, 2007 G3'."
				   "Currently only 2D is supported to calculate the effective strain rate.")
  }
}
