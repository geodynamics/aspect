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



#ifndef __aspect__boundary_temperature_dynamic_core_h
#define __aspect__boundary_temperature_dynamic_core_h

#include <aspect/boundary_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace BoundaryTemperature
  {
    
    /**
     * Data structure for core energy balance calculation
     */
    struct _Core_Data
    {
        /**
         * Energy for specific heat, radioactive heating, gravitational contribution, 
         * adiabatic contribution, and latent heat.
         */
        double Qs,Qr,Qg,Qk,Ql;

        /**
         * Entropy for specific heat, radioactive heating, gravitational contribution, 
         * adiabatic contribution, latent heat, and heat solution.
         */
        double Es,Er,Eg,Ek,El,Eh;

        /**
         * Parameters for core evolution
         * Ri     inner core radius
         * Ti     core-mantle boundary (CMB) temperature
         * Xi     light component concentration in liquid core
         */
        double Ri,Ti,Xi;

        /**
         * Core-mantle boundary heat flux (Q) and core radioactive heating rate (H)
         */
        double Q,H;

        /**
         * Time step for core energy balance solver
         */
        double dt;
        
        /**
         * The changing rate of inner core radius, CMB temperature, and light component
         * concentration.
         */
        double dR_dt,dT_dt,dX_dt;

        /**
         * Other energy source into the core
         */
        double Q_OES;

        bool is_initialized;
    };
    
    
    /**
     * A class that implements a temperature boundary condition for a spherical
     * shell geometry in which the temperature at the outter surfaces are constant
     * and the core-mantle boundaries (CMB) temperature is calculated by core energy balance.
     * The formulation of core energy balance are from: 
     * Nimmo et al., The influence of potassium on core and geodynamo evolution. 
     *    Geophysical Journal International, 2004. 156(2): p. 363-376.
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class Dynamic_core : public Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /** 
         * Construstor
         */
        Dynamic_core();
        
        /**
         * Return the temperature that is to hold at a particular location on the
         * boundary of the domain. This function returns the temperatures
         * at the inner and outer boundaries.
         *
         * @param geometry_model The geometry model that describes the domain. This may
         *   be used to determine whether the boundary temperature model is implemented
         *   for this geometry.
         * @param boundary_indicator The boundary indicator of the part of the boundary
         *   of the domain on which the point is located at which we are requesting the
         *   temperature.
         * @param location The location of the point at which we ask for the temperature.
         **/
        virtual
        double  boundary_temperature (const types::boundary_id            boundary_indicator,
                                      const Point<dim>                    &location) const;

        /**
         * Return the minimal the temperature on that part of the boundary
         * on which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Return the maximal the temperature on that part of the boundary
         * on which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Declare the parameters this class takes through input files.
         * This class declares the inner and outer boundary temperatures.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        
        /**
         * This function update the core-mantle boundary (CMB) temperature by 
         * the core energy balance solver using the core-mantle boundary heat flux.
         */
        virtual void update();

        /**
         * Pass core data to other modules
         */
        struct _Core_Data get_core_data() const;

        /**
         * Check if other energy source in core is in use
         */
        bool is_OES_used() const;


      private:

        /**
         * Data for core energy balance
         * it get updated each time step.
         */
        struct _Core_Data core_data;

        /**
         * Temperatures at the inner boundaries.
         */
        double inner_temperature;
        
        /**
         * Temperatures at the outer boundaries.
         */
        double outer_temperature;

        /**
         * Initial CMB temperature changing rate
         */
        double init_dT_dt;

        /**
         * Initial inner core radius changing rate
         */
        double init_dR_dt;

        /**
         * Initial light composition changing rate
         */
        double init_dX_dt;

        /**
         * Flag for determining the initial call for update().
         */
        bool is_first_call;
        
        /**
         * Gravitation constant
         */
        const double G;
        
        /**
         * Core radius
         */
        double Rc;
        
        /**
         * (Heat capacity) * density
         */
        double CpRho;
        
        /**
         * Initial light composition concentration
         */
        double X_init;

        /**
         * Partition coefficient of the light element
         */
        double Delta;

        /**
         * Gravitational acceleration
         */
        double g;

        /**
         * Pressure at the core mantle boundary
         */
        double P_CMB;

        /**
         * Pressure at the center of the core
         */
        double P_Core;

        /**
         * Parameters for core solidus following:
         * if not dependent on composition 
         *   Tm(p)= Tm0*(1-Theta)*(1+Tm1*p+Tm2*p^2)
         * if depend on composition X
         *   Tm(p)= Tm0*(1-Theta*X)*(1+Tm1*p+Tm2*p^2)
         */
        double Tm0;
        double Tm1;
        double Tm2;
        double Theta;
        bool composition_dependency;
        
        /**
         * If using the Fe-FeS system solidus from Buono & Walker (2011) instead.
         */
        bool use_bw11;

        //Variables for formulation of Nimmo et al. [2004]
        /**
         * Compressibility at zero pressure
         */
        double K0;
        
        /**
         * Thermal expansivity
         */
        double Alpha;
        
        /**
         * Density at zero pressure
         */
        double Rho_0;
        
        /**
         * Density at the center of the planet
         */
        double Rho_cen;
        
        /**
         * Latent heat of fusion
         */
        double Lh;

        /**
         * Heat of reaction
         */
        double Rh;
        
        /**
         * Compositional expansion coefficient
         */
        double Beta_c;
        
        /**
         * Heat conductivity of the core
         */
        double k_c;
        
        /**
         * Heat capacity
         */
        double Cp;

        /**
         * Number of radioheating element in core
         */
        unsigned int n_radioheating_elements;
        
        /**
         * Heating rates of different elements
         */
        std::vector<double> heating_rate;
        
        /**
         * Half life of different elements
         */
        std::vector<double> half_life;
        
        /**
         * Initial concentration of different elements
         */
        std::vector<double> initial_concentration;
        
        /**
         * Two length scales in Nimmo et al. (2004)
         */
        double L;
        double D;
        
        /**
         * Mass of the core
         */
        double Mc;

        /**
         * Max iterations for the core energy balance solver.
         */
        int max_steps;

        /**
         * Temperature correction value for adiabatic
         */
        double dTa;

        /**
         * Other energy source into the core, e.g. the mechanical stirring of the moon.
         */
        std::string name_OES;
        struct str_data_OES
        {
          double t;
          double w;
        };
        std::vector<struct str_data_OES> data_OES;
        void read_data_OES();
        double get_OES(double t) const;

        

        /**
         * Solve core energy balance for each time step.
         * Well solving the change in core-mantle boundary temperature T, inner core radius R, and 
         *    light component (e.g. S, O, Si) composition X, the following relations has to be respected:
         * 1. At the inner core boundary the adiabatic temperature should be equal to solidu temperature
         * 2. The following energy production rate should be balanced in core:
         *    Heat flux at core-mantle boundary         Q
         *    Specific heat                             Qs*dT/dt
         *    Radioactive heating                       Qr
         *    Gravitational contribution                Qg*dR/dt
         *    Latent heat                               Ql*dR/dt
         *    So that         Q+Qs*dT/dt+Qr+Qg*dR/dt*Ql*dR/dt=0
         * 3. The light component composition X depends on inner core radius (See function get_X() ),
         *    and core solidus may dependent on X as well.
         *    This becomes a small nonliner problem. Directly iterate through the above three equations doesn't
         *    converge well. Alternatively we solve the inner core radius by bisection method.
         */
        bool solve_time_step(double &X, double &T, double &R);
        
        /**
         * Compute the difference between solidus and adiabatic at inner 
         * core boundary for a given inner core radius.
         */
        double get_dT(double r) const;

        /**
         * Use energy balance to calculate core mantle boundary temperature 
         * with a given inner core radius.
         */
        double get_Tc(double r) const;

        /**
         * Get the solidus temperature at inner core boundary 
         * with a given inner core radius.
         */
        double get_Ts(double r) const;
        
        /**
         * Compute the core solidus at certain pressure
         */
        double get_solidus(double X,double p) const;

        /**
         * Get initial inner core radius with given initial core mantle temperature. 
         */
        double get_initial_Ri(double T);
        
        /**
         * Get the light composition concentration in the outter core from given
         * inner core radius r
         */
        double get_X(double r) const;
        
        /**
         * Compute the mass inside certain radius with in the core_data.
         */
        double get_Mass(double r) const;
        
        /**
         * Calculate Sn(B,R), referring to Nimmo et al. (2004) paper
         */
        double fun_Sn(double B,double R,double n) const;
        
        /**
         * Calculate density at given r
         */
        double get_Rho(double r) const;
        
        /**
         * Calculate gravitational acceleration at given r
         */
        double get_g(double r) const;
        
        /**
         * Calculate the core temperature at given r
         * Tc is the temperature at CMB
         */
        double get_T(double Tc, double r) const;
        
        /**
         * Calculate pressure at given r
         */
        double get_Pressure(double r) const;
        
        /**
         * Calculate the gravitational potential at given r
         */
        double get_gravity_potential(double r) const;
        
        /**
         * Calculate energy and entropy change rate factor (regarding the core 
         * cooling rated Tc/dt) Qs and Es with given core-mantle boundary (CMB)
         * temperature Tc
         */
        void get_specific_heating(double Tc, double &Qs,double &Es);

        /**
         * Calculate energy and entropy change rate factor (regarding the 
         * radioactive heating rate H) Qr and Er with given CMB temperature Tc
         */
        void get_radio_heating(double Tc, double &Qr, double &Er);
        
        /**
         * Calculate energy and entropy change rate factor (regarding the inner core 
         * growth rate dR/dt) Qg and Eg with given Tc(CMB temperature), r(inner core 
         * radius), X(light composition concentration)
         */
        void get_gravity_heating(double Tc, double r,double X,double &Qg,double &Eg);
        
        /**
         * Calculate energy and entropy change rate factor (regarding the core 
         * cooling rate Tc/dt) Qk and Ek with given Tc(CMB temperature)
         */
        void get_adiabatic_heating(double Tc, double &Ek, double &Qk);
        
        /**
         * Calculate energy and entropy change rate factor (regarding the inner core 
         * growth rate dR/dt) Ql and El with given Tc(CMB temperature), r(inner core 
         * radius)
         */
        void get_latent_heating(double Tc, double r, double &El, double &Ql);

        /**
         * Calculate entropy of heat of solution Eh
         */
        void get_heat_solution(double Tc, double r, double X, double &Eh);
        
        /**
         * return radio heating rate at certain time
         */
        double get_radioheating_rate() const;
        
        /**
         * Update the data for core dynamic simulation, the data will be used   
         * in the next timestep and for postprocess.
         */
        void update_core_data();

    };
  }
}


#endif
