/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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



#include <aspect/boundary_temperature/dynamic_core.h>
#include <aspect/postprocess/core_statistics.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/exceptions.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {

    template <int dim>
    double
    DynamicCore<dim>::
    boundary_temperature (const types::boundary_id            boundary_indicator,
                          const Point<dim>                    &/*location*/) const
    {
      switch (boundary_indicator)
        {
          case 0:
            return inner_temperature;
          case 1:
            return outer_temperature;
          default:
            Assert (false, ExcMessage ("Unknown boundary indicator."));
            return std::numeric_limits<double>::quiet_NaN();
        }
    }


    template <int dim>
    double
    DynamicCore<dim>::
    minimal_temperature (const std::set<types::boundary_id> &/*fixed_boundary_ids*/) const
    {
      return std::min (inner_temperature, outer_temperature);
    }



    template <int dim>
    double
    DynamicCore<dim>::
    maximal_temperature (const std::set<types::boundary_id> &/*fixed_boundary_ids*/) const
    {
      return std::max (inner_temperature, outer_temperature);
    }



    template <int dim>
    void
    DynamicCore<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Dynamic core");
        {
          prm.declare_entry ("Outer temperature", "0.",
                             Patterns::Double (),
                             "Temperature at the outer boundary (lithosphere water/air). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Inner temperature", "6000.",
                             Patterns::Double (),
                             "Temperature at the inner boundary (core mantle boundary) at the "
                             "beginning. Units: \\si{\\kelvin}.");
          prm.declare_entry ("dT over dt", "0.",
                             Patterns::Double (),
                             "Initial CMB temperature changing rate. "
                             "Units: \\si{\\kelvin}/year.");
          prm.declare_entry ("dR over dt", "0.",
                             Patterns::Double (),
                             "Initial inner core radius changing rate. "
                             "Units: \\si{\\kilo\\meter}/year.");
          prm.declare_entry ("dX over dt", "0.",
                             Patterns::Double (),
                             "Initial light composition changing rate. "
                             "Units: 1/year.");
          prm.declare_entry ("Core density", "12.5e3",
                             Patterns::Double (),
                             "Density of the core. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Gravity acceleration", "9.8",
                             Patterns::Double (),
                             "Gravitation acceleration at CMB. "
                             "Units: \\si{\\meter\\per\\second\\squared}.");
          prm.declare_entry ("CMB pressure", "0.14e12",
                             Patterns::Double (),
                             "Pressure at CMB. Units: \\si{\\pascal}.");
          prm.declare_entry ("Initial light composition", "0.01",
                             Patterns::Double (0.),
                             "Initial light composition (eg. S,O) concentration "
                             "in weight fraction.");
          prm.declare_entry ("Max iteration", "30000",
                             Patterns::Integer (0),
                             "The max iterations for nonlinear core energy solver.");
          prm.declare_entry ("Core heat capacity", "840.",
                             Patterns::Double (0.),
                             "Heat capacity of the core. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("K0", "4.111e11",
                             Patterns::Double (0.),
                             "Core compressibility at zero pressure. "
                             "See \\cite{NPB+04} for more details.");
          prm.declare_entry ("Rho0", "7.019e3",
                             Patterns::Double (0.),
                             "Core density at zero pressure. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}. "
                             "See \\cite{NPB+04} for more details.");
          prm.declare_entry ("Alpha", "1.35e-5",
                             Patterns::Double (0.),
                             "Core thermal expansivity. Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Lh", "750e3",
                             Patterns::Double (0.),
                             "The latent heat of core freeze. "
                             "Units: \\si{\\joule\\per\\kilogram}.");
          prm.declare_entry ("Rh","-27.7e6",
                             Patterns::Double (),
                             "The heat of reaction. "
                             "Units: \\si{\\joule\\per\\kilogram}.");
          prm.declare_entry ("Beta composition", "1.1",
                             Patterns::Double (0.),
                             "Compositional expansion coefficient $Beta_c$. "
                             "See \\cite{NPB+04} for more details.");
          prm.declare_entry ("Delta","0.5",
                             Patterns::Double (0., 1.),
                             "Partition coefficient of the light element.");
          prm.declare_entry ("Core conductivity", "60.",
                             Patterns::Double (0.),
                             "Core heat conductivity $k_c$. Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.enter_subsection("Geotherm parameters");
          {
            prm.declare_entry ("Tm0","1695.",
                               Patterns::Double (0.),
                               "Melting curve (\\cite{NPB+04} eq. (40)) parameter Tm0. Units: \\si{\\kelvin}.");
            prm.declare_entry ("Tm1","10.9",
                               Patterns::Double (),
                               "Melting curve (\\cite{NPB+04} eq. (40)) parameter Tm1. "
                               "Units: \\si{\\per\\tera\\pascal}.");
            prm.declare_entry ("Tm2","-8.0",
                               Patterns::Double (),
                               "Melting curve (\\cite{NPB+04} eq. (40)) parameter Tm2. "
                               "Units: \\si{\\per\\tera\\pascal\\squared}.");
            prm.declare_entry ("Theta","0.11",
                               Patterns::Double (),
                               "Melting curve (\\cite{NPB+04} eq. (40)) parameter Theta.");
            prm.declare_entry ("Composition dependency","true",
                               Patterns::Bool (),
                               "If melting curve dependent on composition.");
            prm.declare_entry ("Use BW11","false",
                               Patterns::Bool (),
                               "If using the Fe-FeS system solidus from Buono \\& Walker (2011) instead.");
          }
          prm.leave_subsection ();

          prm.enter_subsection("Radioactive heat source");
          {
            prm.declare_entry ("Number of radioactive heating elements","0",
                               Patterns::Integer (0),
                               "Number of different radioactive heating elements in core");
            prm.declare_entry ("Heating rates","",
                               Patterns::List (Patterns::Double ()),
                               "Heating rates of different elements (W/kg)");
            prm.declare_entry ("Half life times","",
                               Patterns::List (Patterns::Double ()),
                               "Half decay times of different elements (Ga)");
            prm.declare_entry ("Initial concentrations","",
                               Patterns::List (Patterns::Double ()),
                               "Initial concentrations of different elements (ppm)");
          }
          prm.leave_subsection ();

          prm.enter_subsection("Other energy source");
          {
            prm.declare_entry ("File name","",
                               Patterns::Anything(),
                               "Data file name for other energy source into the core. "
                               "The 'other energy source' is used for external core energy source."
                               "For example if someone want to test the early lunar core powered by precession "
                               "(Dwyer, C. A., et al. (2011). A long-lived lunar dynamo driven by continuous mechanical stirring. Nature 479(7372): 212-214.)"
                               "Format [Time(Gyr)   Energy rate(W)]");
          }
          prm.leave_subsection ();

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    DynamicCore<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Dynamic core");
        {
          // verify that the geometry is in fact a spherical shell since only
          // for this geometry do we know for sure what boundary indicators it
          // uses and what they mean
          AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()),
                       ExcMessage ("This boundary model is only implemented if the geometry is "
                                   "a spherical shell."));

          inner_temperature = prm.get_double ("Inner temperature");
          outer_temperature = prm.get_double ("Outer temperature");
          init_dT_dt        = prm.get_double ("dT over dt") / year_in_seconds;
          init_dR_dt        = prm.get_double ("dR over dt") / year_in_seconds * 1.e3;
          init_dX_dt        = prm.get_double ("dX over dt") / year_in_seconds;
          Rho_cen           = prm.get_double ("Core density");
          g                 = prm.get_double ("Gravity acceleration");
          P_CMB             = prm.get_double ("CMB pressure");
          X_init            = prm.get_double ("Initial light composition");
          max_steps         = prm.get_integer ("Max iteration");
          Cp                = prm.get_double ("Core heat capacity");
          CpRho             = Cp*Rho_cen;

          //\cite{NPB+04}
          K0                = prm.get_double ("K0");
          Alpha             = prm.get_double ("Alpha");
          Rho_0             = prm.get_double ("Rho0");
          Lh                = prm.get_double ("Lh");
          Beta_c            = prm.get_double ("Beta composition");
          k_c               = prm.get_double ("Core conductivity");
          Delta             = prm.get_double ("Delta");
          Rh                = prm.get_double ("Rh");

          prm.enter_subsection("Geotherm parameters");
          {
            Tm0           =  prm.get_double ("Tm0");
            Tm1           =  prm.get_double ("Tm1");
            Tm2           =  prm.get_double ("Tm2");
            Theta         =  prm.get_double ("Theta");
            composition_dependency
              =  prm.get_bool("Composition dependency");
            use_bw11      =  prm.get_bool("Use BW11");
          }
          prm.leave_subsection ();

          prm.enter_subsection("Radioactive heat source");
          {
            n_radioheating_elements = prm.get_integer ("Number of radioactive heating elements");
            heating_rate = Utilities::string_to_double
                           (Utilities::split_string_list
                            (prm.get("Heating rates")));
            AssertThrow(n_radioheating_elements==heating_rate.size(),
                        ExcMessage("Number of heating rate entities does not match "
                                   "the number of radioactive elements."));
            half_life = Utilities::string_to_double
                        (Utilities::split_string_list
                         (prm.get("Half life times")));
            AssertThrow(n_radioheating_elements==half_life.size(),
                        ExcMessage("Number of half life time entities does not match "
                                   "the number of radioactive elements."));
            initial_concentration = Utilities::string_to_double
                                    (Utilities::split_string_list
                                     (prm.get("Initial concentrations")));
            AssertThrow(n_radioheating_elements==initial_concentration.size(),
                        ExcMessage("Number of initial concentration entities does not match "
                                   "the number of radioactive elements."));
          }
          prm.leave_subsection ();

          prm.enter_subsection("Other energy source");
          {
            name_OES = prm.get("File name");
          }
          prm.leave_subsection ();

          L=sqrt(3*K0*(log(Rho_cen/Rho_0)+1)/(2*M_PI*constants::big_g*Rho_0*Rho_cen));
          D=sqrt(3*Cp/(2*M_PI*Alpha*Rho_cen*constants::big_g));

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    DynamicCore<dim>::read_data_OES()
    {
      data_OES.clear();
      if (name_OES.size()==0)
        return;
      std::istringstream in(Utilities::read_and_distribute_file_content(name_OES.c_str(),
                                                                        this->get_mpi_communicator()));
      if (in.good())
        {
          str_data_OES data_read;
          std::string line;
          while (!in.eof())
            {
              std::getline(in, line);
              if (sscanf(line.data(), "%le\t%le\n", &data_read.t, &data_read.w)==2)
                data_OES.push_back(data_read);
            }
        }
      if (data_OES.size()!=0)
        this->get_pcout() << "Other energy source is in use ( "
                          << data_OES.size()
                          << " data points is read)."
                          << std::endl;
    }

    template <int dim>
    double
    DynamicCore<dim>::get_OES(double t) const
    {
      // The core evolution is quite slow, so the time units used here is billion years.
      t/=1.e9*year_in_seconds;
      double w=0.;
      for (unsigned i=1; i<data_OES.size(); i++)
        {
          if (t>=data_OES[i-1].t && t<data_OES[i].t )
            {
              w = data_OES[i-1].w + ( t - data_OES[i-1].t)
                  /(data_OES[i].t - data_OES[i-1].t)
                  *(data_OES[i].w - data_OES[i-1].w);
              break;
            }
        }
      return w;
    }

    template <int dim>
    DynamicCore<dim>::DynamicCore()
    {
      is_first_call = true;
      core_data.is_initialized = false;
    }

    template <int dim>
    double
    DynamicCore<dim>::get_initial_Ri(double T)
    {
      double r0=0.,
             r1=Rc;
      double dT0=get_T(T,r0)-get_solidus(get_X(r0),get_Pressure(r0)),
             dT1=get_T(T,r1)-get_solidus(get_X(r1),get_Pressure(r1));
      if (dT0<=0. && dT1<=0.)
        return Rc;
      if (dT0>=0. && dT1>=0.)
        return 0.;
      for (int i=0; i<max_steps; i++)
        {
          double rm=(r0+r1)/2.,
                 dTm=get_T(T,rm)-get_solidus(get_X(rm),get_Pressure(rm));
          if (dTm==0.)
            return rm;
          if (dTm*dT0<0.)
            {
              r1=rm;
              continue;
            }
          if (dTm*dT1<0.)
            {
              r0=rm;
              continue;
            }
        }
      if (dT0>0 && dT1<0)
        {
          // Snowing core
          AssertThrow(false, ExcMessage("[Dynamic core] You had a 'Snowing Core' (i.e., core is freezing from CMB), "
                                        "the treatment is not available at the moment."));
        }
      return (r0+r1)/2.;
    }

    template <int dim>
    bool
    DynamicCore<dim>::solve_time_step(double &X, double &T, double &R)
    {
      // Well solving the change in core-mantle boundary temperature T, inner core radius R, and
      //    light component (e.g. S, O, Si) composition X, the following relations has to be respected:
      // 1. At the inner core boundary the adiabatic temperature should be equal to solidus temperature
      // 2. The following energy production rate should be balanced in core:
      //    Heat flux at core-mantle boundary         Q
      //    Specific heat                             Qs*dT/dt
      //    Radioactive heating                       Qr
      //    Gravitational contribution                Qg*dR/dt
      //    Latent heat                               Ql*dR/dt
      //    So that         Q+Qs*dT/dt+Qr+Qg*dR/dt*Ql*dR/dt=0
      // 3. The light component composition X depends on inner core radius (See function get_X() ),
      //    and core solidus may dependent on X as well
      // This becomes a small nonlinear problem. Directly iterate through the above three system doesn't
      // converge well. Alternatively we solve the inner core radius by bisection method.

      int steps=1;
      double R_0,R_1,R_2;
      // dT is the temperature difference between adiabatic and solidus at
      // inner-outer core boundary. If dT=0 then we found our solution.
      double dT0,dT1,dT2;
      R_0 = 0.;
      R_1 = core_data.Ri;
      R_2 = Rc;
      dT0 = get_dT(R_0);
      dT1 = get_dT(R_1);
      dT2 = get_dT(R_2);

      if (dT0 >= 0. && dT2 >= 0.)
        {
          // Fully molten core
          R_1 = R_0;
          dT1 = 0;
        }
      else if (dT2 <= 0. && dT0 <= 0. )
        {
          //Completely solid core
          R_1 = R_2;
          dT1 = 0;
        }
      else while (!(dT1==0 || steps>max_steps))
          {
            // If solution is out of the interval, then something is wrong.
            if (dT0*dT2>0)
              {
                this->get_pcout()<<"Step: "<<steps<<std::endl
                                 <<" R=["<<R_0/1e3<<","<<R_2/1e3<<"]"<<"(km)"
                                 <<" dT0="<<dT0<<", dT2="<<dT2<<std::endl
                                 <<"Q_CMB="<<core_data.Q<<std::endl
                                 <<"Warning: Solution for inner core radius can not be found! Mid-point is used."<<std::endl;
                AssertThrow(dT0*dT2<=0,ExcMessage("No single solution for inner core!"));
              }
            else if (dT0*dT1<0.)
              {
                R_2 = R_1;
                dT2 = dT1;
              }
            else if (dT2*dT1<0.)
              {
                R_0 = R_1;
                dT0 = dT1;
              }
            R_1 = (R_0 + R_2) / 2.;
            dT1 = get_dT(R_1);
            steps++;
          }

      // Calculate new R,T,X
      R = R_1;
      T = get_Tc(R);
      X = get_X(R);

      if (dT0<0. && dT2>0.)
        {
          // Normal solution
          return true;
        }
      else if (dT0>0. && dT2<0.)
        {
          // Snowing core solution
          return false;
        }
      else
        {
          // No solution found.
          this->get_pcout()<<"[Dynamic core] Step: "<<steps<<std::endl
                           <<" R=["<<R_0/1e3<<","<<R_2/1e3<<"]"<<"(km)"
                           <<" dT0="<<dT0<<", dT2="<<dT2<<std::endl
                           <<"Q_CMB="<<core_data.Q<<std::endl;
          AssertThrow(false, ExcMessage("[Dynamic core] No inner core radius solution found!"));
        }

      return false;
    }

    template <int dim>
    double
    DynamicCore<dim>::get_Tc(double r) const
    {
      // Using all Q values from last step.
      // Qs & Qr is constant, while Qg & Ql depends on inner core radius Ri
      // TODO: Use mid-point value for Q values.
      return core_data.Ti - ( (core_data.Q + core_data.Qr + core_data.Q_OES) * core_data.dt
                              + (core_data.Qg + core_data.Ql)*(r-core_data.Ri)
                            ) / core_data.Qs;
    }

    template <int dim>
    double
    DynamicCore<dim>::get_Ts(double r) const
    {
      return get_solidus(get_X(r),get_Pressure(r));
    }


    template <int dim>
    double
    DynamicCore<dim>::get_dT(double r) const
    {
      return get_T(get_Tc(r),r) - get_Ts(r);
    }

    template <int dim>
    void
    DynamicCore<dim>::update_core_data()
    {
      get_specific_heating(core_data.Ti,core_data.Qs,core_data.Es);
      get_radio_heating(core_data.Ti,core_data.Qr,core_data.Er);
      get_gravity_heating(core_data.Ti,core_data.Ri,core_data.Xi,core_data.Qg,core_data.Eg);
      get_adiabatic_heating(core_data.Ti,core_data.Ek,core_data.Qk);
      get_latent_heating(core_data.Ti,core_data.Ri,core_data.El,core_data.Ql);
      get_heat_solution(core_data.Ti,core_data.Ri,core_data.Xi,core_data.Eh);
    }

    template <int dim>
    const internal::CoreData &
    DynamicCore<dim>::get_core_data() const
    {
      return core_data;
    }

    template <int dim>
    double
    DynamicCore<dim>::get_solidus(double X,double p) const
    {
      if (use_bw11)
        {
          // Change from weight percent to mole percent.
          double x,x0=32./88.;
          if (X<x0)
            x=56.*X/(32.*(1.-X));
          else
            x=1.;
          // Change from Pa to GPa
          p*=1e-9;
          // Fe-FeS system solidus by Buono & Walker (2011)
          return (-2.4724*p*p*p*p + 28.025*p*p*p + 9.1404*p*p + 581.71*p + 3394.8) * x*x*x*x
                 +( 1.7978*p*p*p*p - 6.7881*p*p*p - 197.69*p*p - 271.69*p - 8219.5) * x*x*x
                 +(-0.1702*p*p*p*p - 9.3959*p*p*p + 163.53*p*p - 319.35*p + 5698.6) * x*x
                 +(-0.2308*p*p*p*p + 7.1000*p*p*p - 64.118*p*p + 105.98*p - 1621.9) * x
                 +( 0.2302*p*p*p*p - 5.3688*p*p*p + 38.124*p*p - 46.681*p + 1813.8);

        }
      else
        {
          if (composition_dependency)
            return (Tm0*(1-Theta*X)*(1+Tm1*p+Tm2*pow(p,2)));
          else
            return (Tm0*(1-Theta)*(1+Tm1*p+Tm2*pow(p,2)));
        }
    }

    template <int dim>
    double
    DynamicCore<dim>::get_X(double r) const
    {
      double xi_3=pow(r/Rc,3);
      return X_init/(1-xi_3+Delta*xi_3);
    }

    template <int dim>
    void
    DynamicCore<dim>::update()
    {
      core_data.dt    = this->get_timestep();
      core_data.H     = get_radioheating_rate();

      // It's a bit tricky here.
      // Didn't use the initialize() function instead because the postprocess is initialized after boundary temperature.
      // It is not available at the time initialize() function of boundary temperature is called.
      if (is_first_call==true)
        {
          AssertThrow(this->get_postprocess_manager().template has_matching_postprocessor<const Postprocess::CoreStatistics<dim> >(),
                      ExcMessage ("Dynamic core boundary condition has to work with dynamic core statistics postprocessor."));

          const Postprocess::CoreStatistics<dim> &core_statistics
            = this->get_postprocess_manager().template get_matching_postprocessor<const Postprocess::CoreStatistics<dim> >();
          // The restart data is stored in 'core statistics' postprocessor.
          // If restart from checkpoint, extract data from there.
          core_data = core_statistics.get_core_data();

          // Read data of other energy source
          read_data_OES();

          const GeometryModel::SphericalShell<dim> &spherical_shell_geometry =
            Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());

          Rc=spherical_shell_geometry.inner_radius();
          Mc=get_Mass(Rc);
          P_Core=get_Pressure(0);

          // If the material model is incompressible, we have to get correction for the real core temperature
          if (this->get_adiabatic_conditions().is_initialized() && !this->get_material_model().is_compressible())
            {
              Point<dim> p1;
              p1(0) = spherical_shell_geometry.inner_radius();
              dTa   = this->get_adiabatic_conditions().temperature(p1)
                      - this->get_adiabatic_surface_temperature();
            }
          else
            dTa   = 0.;

          // Setup initial core data from prm input.
          // If resumed from checkpoint, core_data is read from postprocess instead of set from prm file.
          // (The boundary_temperature doesn't seem to support restart/resume, the data has to passed and
          // stored in the postprocessor 'core statistics')
          if (!core_data.is_initialized)
            {
              core_data.Ti=inner_temperature + dTa;
              core_data.Ri = get_initial_Ri(core_data.Ti);
              core_data.Xi = get_X(core_data.Ri);

              core_data.Q=0.;
              core_data.dt=0.;
              core_data.dT_dt=init_dT_dt;
              core_data.dR_dt=init_dR_dt;
              core_data.dX_dt=init_dX_dt;
              update_core_data();
              core_data.is_initialized = true;
              std::stringstream output;
              output<<std::setiosflags(std::ios::left)
                    <<"   Dynamic core initialized as:"<<std::endl
                    <<"     "<<std::setw(15)<<"Tc(K)"<<std::setw(15)<<"Ri(km)"<<std::setw(15)<<"Xi"
                    <<std::setw(15)<<"dT/dt(K/year)"<<std::setw(15)<<"dR/dt(km/year)"<<std::setw(15)<<"dX/dt(1/year)"<<std::endl
                    <<"     "<<std::setprecision(6)<<std::setw(15)<<inner_temperature<<std::setw(15)<<core_data.Ri/1.e3<<std::setw(15)<<core_data.Xi
                    <<std::setw(15)<<core_data.dT_dt *year_in_seconds<<std::setw(15)<<core_data.dR_dt/1.e3 *year_in_seconds
                    <<std::setw(15)<<core_data.dX_dt *year_in_seconds<<std::endl;
              this->get_pcout() << output.str();
            }
          is_first_call = false;
        }

      // Calculate core mantle boundary heat flow
      {
        const QGauss<dim-1> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
        FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                          this->get_fe(),
                                          quadrature_formula,
                                          update_gradients      | update_values |
                                          update_normal_vectors |
                                          update_quadrature_points       | update_JxW_values);

        std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());
        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        //std::map<types::boundary_id, double> local_boundary_fluxes;
        double local_CMB_flux   = 0.;
        double local_CMB_area   = 0.;

        types::boundary_id CMB_id = 0;

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

        // for every surface face on which it makes sense to compute a
        // heat flux and that is owned by this processor,
        // integrate the normal heat flux given by the formula
        //   j =  - k * n . grad T
        //
        // for the spherical shell geometry, note that for the inner boundary,
        // the normal vector points *into* the core, i.e. we compute the flux
        // *out* of the mantle, not into it. we fix this when we add the local
        // contribution to the global flux

        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->at_boundary(f))
                if (cell->face(f)->boundary_id() == CMB_id)
                  {
                    fe_face_values.reinit (cell, f);
                    fe_face_values[this->introspection().extractors.temperature].get_function_gradients (this->get_solution(),
                        temperature_gradients);
                    fe_face_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                        in.temperature);
                    fe_face_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                                   in.pressure);
                    for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      fe_face_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                          composition_values[c]);

                    in.position = fe_face_values.get_quadrature_points();

                    // since we are not reading the viscosity and the viscosity
                    // is the only coefficient that depends on the strain rate,
                    // we need not compute the strain rate. set the corresponding
                    // array to empty, to prevent accidental use and skip the
                    // evaluation of the strain rate in evaluate().
                    in.strain_rate.resize(0);

                    for (unsigned int i=0; i<fe_face_values.n_quadrature_points; ++i)
                      {
                        for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                          in.composition[i][c] = composition_values[c][i];
                      }

                    this->get_material_model().evaluate(in, out);


                    double local_normal_flux = 0;
                    double local_face_area   = 0;
                    for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                      {
                        const double thermal_conductivity
                          = out.thermal_conductivities[q];
                        double adiabatic_flux=0.;
                        if (this->get_material_model().is_compressible()==false)
                          {
                            const double alpha = out.thermal_expansion_coefficients[q];
                            const double cp = out.specific_heat[0];
                            const double gravity = this->get_gravity_model().gravity_vector(in.position[q]).norm();
                            if (cell->face(f)->boundary_id()==0)
                              adiabatic_flux = - alpha * gravity / cp;
                            else if (cell->face(f)->boundary_id()==1)
                              adiabatic_flux = alpha * gravity / cp;
                          }

                        local_normal_flux += -thermal_conductivity *
                                             (temperature_gradients[q] * fe_face_values.normal_vector(q)
                                              + adiabatic_flux) * fe_face_values.JxW(q);
                        local_face_area   += fe_face_values.JxW(q);

                      }
                    local_CMB_flux += local_normal_flux;
                    local_CMB_area += local_face_area;
                  }
        // now communicate to get the global values
        double global_CMB_flux;
        double global_CMB_area;
        global_CMB_flux = Utilities::MPI::sum (local_CMB_flux, this->get_mpi_communicator());
        global_CMB_area = Utilities::MPI::sum (local_CMB_area, this->get_mpi_communicator());

        // Using area averaged heat-flux density times core mantle boundary area to calculate total heat-flux on the 3D sphere.
        // By doing this, using dynamic core evolution with geometry other than 3D spherical shell becomes possible.
        double average_CMB_heatflux_density = global_CMB_flux / global_CMB_area;
        core_data.Q = average_CMB_heatflux_density * 4. * M_PI * Rc * Rc;
      }

      core_data.Q_OES = get_OES(this->get_time());

      if ((core_data.Q + core_data.Q_OES) * core_data.dt!=0.)
        {
          double X1,R1=core_data.Ri,T1;
          solve_time_step(X1,T1,R1);
          if (core_data.dt!=0)
            {
              core_data.dR_dt=(R1-core_data.Ri)/core_data.dt;
              core_data.dT_dt=(T1-core_data.Ti)/core_data.dt;
              core_data.dX_dt=(X1-core_data.Xi)/core_data.dt;
            }
          else
            {
              core_data.dR_dt=0.;
              core_data.dT_dt=0.;
              core_data.dX_dt=0.;
            }
          core_data.Xi=X1;
          core_data.Ri=R1;
          core_data.Ti=T1;
        }

      inner_temperature = core_data.Ti - dTa;
      update_core_data();
      if ((core_data.Q + core_data.Q_OES + core_data.Qr) * core_data.dt!=0.)
        {
          std::stringstream output;
          output<<std::setiosflags(std::ios::left)
                <<"   Dynamic core data updated."<<std::endl
                <<"     "<<std::setw(15)<<"Tc(K)"<<std::setw(15)<<"Ri(km)"<<std::setw(15)<<"Xi"
                <<std::setw(15)<<"dT/dt(K/year)"<<std::setw(15)<<"dR/dt(km/year)"<<std::setw(15)<<"dX/dt(1/year)"<<std::endl
                <<"     "<<std::setprecision(6)<<std::setw(15)<<inner_temperature<<std::setw(15)<<core_data.Ri/1.e3<<std::setw(15)<<core_data.Xi
                <<std::setw(15)<<core_data.dT_dt *year_in_seconds<<std::setw(15)<<core_data.dR_dt/1.e3 *year_in_seconds
                <<std::setw(15)<<core_data.dX_dt *year_in_seconds<<std::endl;
          this->get_pcout() << output.str();
        }
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_Mass(double r) const
    {
      return 4.*M_PI*Rho_cen*(-pow(L,2)/2.*r*exp(-pow(r/L,2))+pow(L,3)/4.*sqrt(M_PI)*erf(r/L));
    }

    template <int dim>
    double
    DynamicCore<dim>::
    fun_Sn(double B,double R,double n) const
    {
      double S=R/(2.*sqrt(M_PI));
      for (unsigned i=1; i<=n; i++)
        S+=B/sqrt(M_PI)*exp(-pow(i,2)/4.)/i*sinh(i*R/B);
      return S;
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_Pressure(double r) const
    {
      return P_CMB-(4*M_PI*constants::big_g*pow(Rho_cen,2))/3
             *((3*pow(r,2)/10.-pow(L,2)/5)*exp(-pow(r/L,2))
               -(3*pow(Rc,2)/10-pow(L,2)/5)*exp(-pow(Rc/L,2)));
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_Rho(double r) const
    {
      return Rho_cen*exp(-pow(r/L,2));
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_g(double r) const
    {
      return (4*M_PI/3)*constants::big_g*Rho_cen*r*(1-3*pow(r,2)/(5*pow(L,2)));
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_T(double Tc,double r) const
    {
      return Tc*exp((pow(Rc,2)-pow(r,2))/pow(D,2));
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_gravity_potential(double r) const
    {
      return 2./3.*M_PI*constants::big_g*Rho_cen*(pow(r,2)*(1.-3.*pow(r,2)
                                                            /(10.*pow(L,2)))-pow(Rc,2)*(1.-3.*pow(Rc,2)/(10.*pow(L,2))));
    }

    template <int dim>
    void
    DynamicCore<dim>::
    get_specific_heating(double Tc, double &Qs,double &Es)
    {
      double A=sqrt(1./(pow(L,-2)+pow(D,-2)));
      double Is=4.*M_PI*get_T(Tc,0.)*Rho_cen*(-pow(A,2)*Rc/2.*exp(-pow(Rc/A,2))+pow(A,3)*sqrt(M_PI)/4.*erf(Rc/A));

      Qs=-Cp/Tc*Is;
      Es=Cp/Tc*(Mc-Is/Tc);
    }

    template <int dim>
    void
    DynamicCore<dim>::
    get_radio_heating(double Tc, double &Qr, double &Er)
    {
      double B,It;
      if (D>L)
        {
          B=sqrt(1/(1/pow(L,2)-1/pow(D,2)));
          It=4*M_PI*Rho_cen/get_T(Tc,0)*(-pow(B,2)*Rc/2*exp(-pow(Rc/B,2))+pow(B,3)/sqrt(M_PI)/4*erf(Rc/B));
        }
      else
        {
          B=sqrt(1/(pow(D,-2)-pow(L,-2)));
          It=4*M_PI*Rho_cen/get_T(Tc,0)*(pow(B,2)*Rc/2*exp(pow(Rc/B,2))-pow(B,2)*fun_Sn(B,Rc,100)/2);
        }
      Qr=Mc*core_data.H;
      Er=(Mc/Tc-It)*core_data.H;

    }

    template <int dim>
    void
    DynamicCore<dim>::
    get_heat_solution(double Tc, double r, double X,double &Eh)
    {
      double B,It;
      if (D>L)
        {
          B=sqrt(1/(1/pow(L,2)-1/pow(D,2)));
          It=4*M_PI*Rho_cen/get_T(Tc,0)*(-pow(B,2)*Rc/2*exp(-pow(Rc/B,2))+pow(B,3)/sqrt(M_PI)/4*erf(Rc/B));
          It-=4*M_PI*Rho_cen/get_T(Tc,0)*(-pow(B,2)*r/2*exp(-pow(r/B,2))+pow(B,3)/sqrt(M_PI)/4*erf(r/B));
        }
      else
        {
          B=sqrt(1/(pow(D,-2)-pow(L,-2)));
          It=4*M_PI*Rho_cen/get_T(Tc,0)*(pow(B,2)*Rc/2*exp(pow(Rc/B,2))-pow(B,2)*fun_Sn(B,Rc,100)/2);
          It-=4*M_PI*Rho_cen/get_T(Tc,0)*(pow(B,2)*r/2*exp(pow(r/B,2))-pow(B,2)*fun_Sn(B,r,100)/2);
        }
      double Cc=4*M_PI*pow(r,2)*get_Rho(r)*X/(Mc-get_Mass(r));
      Eh=Rh*(It-(Mc-get_Mass(r))/get_T(Tc,r))*Cc;
    }

    template <int dim>
    void
    DynamicCore<dim>::
    get_gravity_heating(double Tc, double r,double X,double &Qg,double &Eg)
    {
      double Cc=4*M_PI*pow(r,2)*get_Rho(r)*X/(Mc-get_Mass(r));
      double C_2=3./16.*pow(L,2)-0.5*pow(Rc,2)*(1.-3./10.*pow(Rc/L,2));
      if (r==Rc)
        Qg=0.;
      else
        Qg=(8./3.*pow(M_PI*Rho_cen,2)*constants::big_g*(
              ((3./20.*pow(Rc,5)-pow(L,2)*pow(Rc,3)/8.-C_2*pow(L,2)*Rc)*exp(-pow(Rc/L,2))
               +C_2/2.*pow(L,3)*sqrt(M_PI)*erf(Rc/L))
              -((3./20.*pow(r,5)-pow(L,2)*pow(r,3)/8.-C_2*pow(L,2)*r)*exp(-pow(r/L,2))
                +C_2/2.*pow(L,3)*sqrt(M_PI)*erf(r/L)))
            -(Mc-get_Mass(r))*get_gravity_potential(r))*Beta_c*Cc;
      Eg=Qg/Tc;

    }

    template <int dim>
    void
    DynamicCore<dim>::
    get_adiabatic_heating(double Tc, double &Ek, double &Qk)
    {
      Ek=16*M_PI*k_c*pow(Rc,5)/5/pow(D,4);
      Qk=8*M_PI*pow(Rc,3)*k_c*Tc/pow(D,2);
    }
    template <int dim>
    void
    DynamicCore<dim>::
    get_latent_heating(double Tc, double r, double &El, double &Ql)
    {
      Ql=4.*M_PI*pow(r,2)*Lh*get_Rho(r);
      El=Ql*(get_T(Tc,r)-Tc)/(Tc*get_T(Tc,r));
    }

    template <int dim>
    double
    DynamicCore<dim>::
    get_radioheating_rate() const
    {
      double time=this->get_time()+0.5*this->get_timestep();
      double Ht=0;
      for (unsigned i=0; i<n_radioheating_elements; i++)
        Ht+=heating_rate[i]*initial_concentration[i]*1e-6*pow(0.5,time/half_life[i]/year_in_seconds/1e9);
      return Ht;
    }

    template <int dim>
    bool
    DynamicCore<dim>::
    is_OES_used() const
    {
      if (data_OES.size()>0)
        return true;
      else
        return false;
    }
  }

}


// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(DynamicCore,
                                               "dynamic core",
                                               "This is a boundary temperature model working only with spherical "
                                               "shell geometry and core statistics postprocessor. The temperature "
                                               "at the top is constant, and the core mantle boundary temperature "
                                               "is dynamically evolving through time by calculating the heat flux "
                                               "into the core and solving the core energy balance. The formulation "
                                               "is mainly following \\cite{NPB+04}, and the plugin is used in "
                                               "Zhang et al. [2016]. The energy of core cooling and freeing of the "
                                               "inner core is included in the plugin. However, current plugin can not "
                                               "deal with the energy balance if the core is in the `snowing core' regime "
                                               "(i.e., the core solidifies from the top instead of bottom)."
                                              )
  }
}
