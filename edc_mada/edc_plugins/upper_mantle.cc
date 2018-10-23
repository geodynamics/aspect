/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


/*
 * include "edge_driven.h"
 */

#include "upper_mantle.h"
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_temperature/adiabatic_boundary.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>


using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    std::vector<double>
    UpperMantle<dim>::
    compute_volume_fractions( const std::vector<double> &compositional_fields) const
    {
          std::vector<double> volume_fractions( compositional_fields.size()+1);

          // clip the compositional fields so they are between zero and one
          std::vector<double> x_comp = compositional_fields;
          for ( unsigned int i=0; i < x_comp.size(); ++i)
            x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);

          // sum the compositional fields for normalization purposes
          double sum_composition = 0.0;
          for ( unsigned int i=0; i < x_comp.size(); ++i)
            sum_composition += x_comp[i];

          if (sum_composition >= 1.0)
            {
              volume_fractions[0] = 0.0;  // background mantle
              for ( unsigned int i=1; i <= x_comp.size(); ++i)
                volume_fractions[i] = x_comp[i-1]/sum_composition;
            }
          else
            {
              volume_fractions[0] = 1.0 - sum_composition; // background mantle
              for ( unsigned int i=1; i <= x_comp.size(); ++i)
                volume_fractions[i] = x_comp[i-1];
            }
          return volume_fractions;
     }

    template <int dim>
    double
    UpperMantle<dim>::
    average_value ( const std::vector<double> &composition,
                    const std::vector<double> &parameter_values,
                    const enum averaging_scheme &average_type) const
    {
      double averaged_parameter = 0.0;
      const std::vector<double> volume_fractions = compute_volume_fractions(composition);

      switch (average_type)
        {
          case arithmetic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*parameter_values[i];
            break;
          }
          case harmonic:
          {
            for (unsigned int i=0; i< volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]/(parameter_values[i]);
            averaged_parameter = 1.0/averaged_parameter;
            break;
          }
          case geometric:
          {
            for (unsigned int i=0; i < volume_fractions.size(); ++i)
              averaged_parameter += volume_fractions[i]*std::log(parameter_values[i]);
            averaged_parameter = std::exp(averaged_parameter);
            break;
          }
          case maximum_composition:
          {
            const unsigned int i = (unsigned int)(std::max_element( volume_fractions.begin(),
                                                                    volume_fractions.end() )
                                                  - volume_fractions.begin());
            averaged_parameter = parameter_values[i];
            break;
          }
          default:
          {
            AssertThrow( false, ExcNotImplemented() );
            break;
          }
        }
      return averaged_parameter;
    }
    template <int dim>
    double
    UpperMantle<dim>::
    diffusion_creep (const double &pressure,
                     const double &temperature) const
    {
    	return  0.5 * std::pow(prefactor_diffusion,-1/stress_exponent_diffusion) *
    	              std::exp((activation_energie_diffusion + pressure*activation_volume_diffusion)/
    	              (constants::gas_constant*temperature*stress_exponent_diffusion)) *
    	              std::pow(grain_size, grain_size_exponent_diffusion);
    }

    template <int dim>
    double
    UpperMantle<dim>::
    dislocation_creep (const double &pressure,
                       const double &temperature,
                       const SymmetricTensor<2,dim> &strain_rate) const
    {

    	const double edot_ii = ( (this->get_timestep_number() == 0 && strain_rate.norm() <= std::numeric_limits<double>::min())
    	                               ?
    	                               ref_strain_rate
    	                               :
    	                               std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate)))),
    	                                        min_strain_rate) );

        return 0.5 * std::pow(prefactor_dislocation,-1/stress_exponent_dislocation) *
                     std::exp((activation_energie_dislocation + pressure*activation_volume_dislocation)/
                     (constants::gas_constant*temperature*stress_exponent_dislocation)) *
                     std::pow(edot_ii,((1. - stress_exponent_dislocation)/stress_exponent_dislocation));
    }
    
    template <>
    void
    UpperMantle<2>::evaluate(const MaterialModel::MaterialModelInputs<2> &in,
                            MaterialModel::MaterialModelOutputs<2> &out) const
    {
      Assert (false, ExcNotImplemented());
      //return 0;
    }

    template <int dim>
    void
    UpperMantle<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // evaluate the base model to get the crustal viscosity
      MaterialModel::MaterialModelOutputs<dim> crust_out = out;
      base_model->evaluate(in, crust_out);

      for (unsigned int i=0; i < in.position.size(); ++i)
        {
          const double temperature = in.temperature[i];
          const double pressure= in.pressure[i];
          const Point<3> pos = in.position[i];
          const std::vector<double> &composition = in.composition[i];
          const double c = (in.composition[i].size()>0)
                           ?
                           std::max(0.0, in.composition[i][0])
                           :
                           0.0;
          // The viscosity of the crust is determined by the base model.
          // The rheology of the lithosphere is governed by disolcation creep flow law
          // and diffusion creep for of the sublithospheric mantle.
		  std::vector<double> composition_viscosities (composition.size()+1);
		  std::vector<double> composition_densities (composition.size()+1);

		  AssertThrow(this->introspection().compositional_name_exists("crust"),
		              ExcMessage("Material model Upper mantle works only works if there is a "
		                         "compositional field called crust."));

		  AssertThrow(this->introspection().compositional_name_exists("mantle_lithosphere"),
		 		      ExcMessage("Material model Upper mantle works only works if there is a "
		 		                  "compositional field called mantle_lithosphere."));
		  const unsigned int crust_idx = this->introspection().compositional_index_for_name("crust")+1;
		  const unsigned int mantle_lithosphere_idx = this->introspection().compositional_index_for_name("mantle_lithosphere")+1;

          const std::vector<double> volume_fractions = compute_volume_fractions(composition);
          std::vector<double> ref_densities(volume_fractions.size(), 3300.0);
          double density = 0.0;

          for (unsigned int c=0; c <= in.composition[i].size() ; ++c)
          {
            if (c == crust_idx)
            {
            	ref_densities[crust_idx] = 2700.0;
            	composition_viscosities[c] =  crust_out.viscosities[i];
            }
            else if (c == mantle_lithosphere_idx)
            {
            	ref_densities[c] = 3300.0;
            	composition_viscosities[c] = std::min(std::max(dislocation_creep (pressure, temperature, in.strain_rate[i]), min_visc), max_visc);
            }
            else
            {
            	ref_densities[c] = 3300.0;
               	double disl  = std::min(std::max(dislocation_creep (pressure, temperature, in.strain_rate[i]), min_visc), max_visc);      
                double diff = std::min(std::max(diffusion_creep (pressure, temperature), min_visc), max_visc);
            	composition_viscosities[c] = (disl * diff) / (disl + diff);
            }
           }

          const double temperature_factor = (1 - thermal_alpha * (in.temperature[i] - reference_T));
                for (unsigned int j=0; j < volume_fractions.size(); ++j)
                     density += volume_fractions[j] * ref_densities[j] * temperature_factor;
          out.densities[i] = density;
          out.viscosities[i] = average_value(composition, composition_viscosities, viscosity_averaging);

          out.thermal_expansion_coefficients[i] = thermal_alpha;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = k_value;
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
    double
    UpperMantle<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    double
    UpperMantle<dim>::
    reference_density () const
    {
      return reference_rho;
    }

    template <int dim>
    double
    UpperMantle<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return thermal_alpha;
    }

    template <int dim>
    double
    UpperMantle<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    double
    UpperMantle<dim>::
    reference_thermal_diffusivity () const
    {
      return k_value/(reference_rho*reference_specific_heat);
    }

    template <int dim>
    bool
    UpperMantle<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    UpperMantle<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Upper mantle");
        {
          prm.declare_entry ("Reference strain rate","1.0e-15",Patterns::Double(0),
                              "Reference strain rate for first time step. Units: $1 / s$");
          prm.declare_entry ("Reference density", "3400",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Reference temperature", "293",
                             Patterns::Double (0),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: $K$.");
          prm.declare_entry ("Reference viscosity", "1e22", Patterns::Double(0),
                             "Reference viscosity for nondimensionalization. Units $Pa s$");
          prm.declare_entry ("Base model","simpler",
                             Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                             "The name of a material model that will be used for the crustal viscosity. "
                             "Valid values for this parameter "
                             "are the names of models that are also valid for the "
                             "``Material models/Model name'' parameter. See the documentation for "
                             "that for more information.");
          prm.declare_entry ("Composition viscosity prefactor", "1.0",
                             Patterns::Double (0),
                             "A linear dependency of viscosity on the first compositional field. "
                             "Dimensionless prefactor. With a value of 1.0 (the default) the "
                             "viscosity does not depend on the composition. See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\xi$ there.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
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
                             "The value of the thermal expansion coefficient $\\alpha$. "
                             "Units: $1/K$.");
          prm.declare_entry ("Density differential for compositional field 1", "0",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
                             "change in composition}$.");
          prm.declare_entry ("Denisity adiabatic condition","true",
                             Patterns::Bool (),
                             "Use adiabatic condition for density");

          //rheology parameters
          prm.declare_entry ("Grain size", "1e-3", Patterns::Double(0), "Units: $m$");
          prm.declare_entry ("Minimum strain rate", "1.4e-18", Patterns::Double(0),
                             "Stabilizes strain dependent viscosity. Units: $1 / s$");
          prm.declare_entry ("Maximum viscosity", "1e34", Patterns::Double(0),
                             "Upper cutoff for effective viscosity. Units: $Pa s$");
          prm.declare_entry ("Minimum viscosity", "1e17", Patterns::Double(0),
                             "Lower cutoff for effective viscosity. Units: $Pa s$");

          //difusion creep parameters
          prm.declare_entry ("Activation energie for diffusion creep", "300e3",
                             Patterns::List(Patterns::Double(0)),
                             "Aactivation energies, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Activation volume for diffusion creep", "6e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponent for diffusion creep", "1",
                             Patterns::List(Patterns::Double(0)),
                             "$n_dislocation$, for background mantle and compositional fields, "
                             "Units: None");
          prm.declare_entry ("Grain size exponent for diffusion creep", "2.5",
                             Patterns::List(Patterns::Double(0)),
                             "$m_diffusion$, for background mantle "
                             "Units: None");
          prm.declare_entry ("Prefactor for diffusion creep", "0.11e-15",
                             Patterns::List(Patterns::Double(0)),
                             "Viscosity prefactor, $A$, for background mantle,  "
                             "Units: $Pa^{-n_{diffusion}} m^{n_{diffusion}/m_{diffusion}} s^{-1}$");
          
          //dislocation creep parameters
          prm.declare_entry ("Activation energie for dislocation creep", "430e3",
                             Patterns::List(Patterns::Double(0)),
                             "Aactivation energies, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Activation volume for dislocation creep", "10e-6",
                             Patterns::List(Patterns::Double(0)),
                             "List of activation volumes, $V_a$, for background mantle Units: $m^3 / mol$");
          prm.declare_entry ("Stress exponent for dislocation creep", "3",
                             Patterns::List(Patterns::Double(0)),
                             "Stress exponent, $n_dislocation$, for background mantle, "
                             "Units: None");
          prm.declare_entry ("Prefactor for dislocation creep", "1.1e-16",
                             Patterns::List(Patterns::Double(0)),
                             "Viscosity prefactor, $A$, for background mantle, "
                             "Units: $Pa^{-n_{dislocation}} m^{n_{dislocation}/m_{dislocation}} s^{-1}$");
          prm.declare_entry ("Water content", "1000",
                             Patterns::List(Patterns::Double(0)),
                             "Water conten for dry olivine, $E_a$, for background mantle Units: $J / mol$");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "When more than one compositional field is present at a point "
                             "with different viscosities, we need to come up with an average "
                             "viscosity at that point.  Select a weighted harmonic, arithmetic, "
                             "geometric, or maximum composition.");
          prm.declare_entry ("Reference compressibility", "4e-12",
                             Patterns::Double (0),
                             "The value of the reference compressibility. "
                             "Units: $1/Pa$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    UpperMantle<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Upper mantle");
        {
          reference_rho                   = prm.get_double ("Reference density");
          reference_T                     = prm.get_double ("Reference temperature");
          eta                             = prm.get_double ("Reference viscosity");

          // Crustal viscosity material model
          AssertThrow( prm.get("Base model") != "edge driven",
                      ExcMessage("You may not use ``edge driven'' as the base model for "
                                  "a this material model.") );
          // create the crustal viscosity model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
             sim->initialize_simulator (this->get_simulator());

          composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
          thermal_viscosity_exponent      = prm.get_double ("Thermal viscosity exponent");
          k_value                         = prm.get_double ("Thermal conductivity");
          reference_specific_heat         = prm.get_double ("Reference specific heat");
          thermal_alpha                   = prm.get_double ("Thermal expansion coefficient");
          compositional_delta_rho         = prm.get_double ("Density differential for compositional field 1");

          //rheology parameters
          grain_size                      = prm.get_double("Grain size");
          min_strain_rate                 = prm.get_double("Minimum strain rate");
          ref_strain_rate                 = prm.get_double("Reference strain rate");
          max_visc                        = prm.get_double ("Maximum viscosity");
          min_visc                        = prm.get_double ("Minimum viscosity");

          //diffusion creep parameters
          activation_energie_diffusion    = prm.get_double ("Activation energie for diffusion creep");
          activation_volume_diffusion     = prm.get_double ("Activation volume for diffusion creep");
          stress_exponent_diffusion       = prm.get_double ("Stress exponent for diffusion creep");
          grain_size_exponent_diffusion   = prm.get_double ("Grain size exponent for diffusion creep");
          prefactor_diffusion             = prm.get_double ("Prefactor for diffusion creep");
          C_OH                            = prm.get_double ("Water content");
          
          //diffusion creep parameters;
          activation_energie_dislocation  = prm.get_double ("Activation energie for dislocation creep");
          activation_volume_dislocation   = prm.get_double ("Activation volume for dislocation creep");
          stress_exponent_dislocation     = prm.get_double ("Stress exponent for dislocation creep");
          prefactor_dislocation           = prm.get_double ("Prefactor for dislocation creep");
          C_OH                            = prm.get_double ("Water content");
          reference_compressibility  = prm.get_double ("Reference compressibility");

          // Rheological parameters
          if (prm.get ("Viscosity averaging scheme") == "harmonic")
              viscosity_averaging = harmonic;
          else if (prm.get ("Viscosity averaging scheme") == "arithmetic")
              viscosity_averaging = arithmetic;
          else if (prm.get ("Viscosity averaging scheme") == "geometric")
              viscosity_averaging = geometric;
          else if (prm.get ("Viscosity averaging scheme") == "maximum composition")
              viscosity_averaging = maximum_composition;
          else
              AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

          density_averaging = arithmetic;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::none;

      this->model_dependence.viscosity |= NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::strain_rate;
      if (composition_viscosity_prefactor != 1.0)
        this->model_dependence.viscosity |= NonlinearDependence::compositional_fields;

      if (thermal_alpha != 0)
        this->model_dependence.density |=NonlinearDependence::temperature;
      if (compositional_delta_rho != 0)
        this->model_dependence.density |=NonlinearDependence::compositional_fields;

      // Parse the parameters of the base model
      base_model->parse_parameters(prm);
      this->model_dependence.viscosity = base_model->get_model_dependence().viscosity;
    }
  }
}
// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(UpperMantle,
                                   "upper mantle",
                                   "A material model for edge driven convection model. There reference density, "
                                   "of the crust is 2700 $kg/m^3$ and 3400 $kg/m^3$ for the mantle. "
                                   "The rheology of the lithosphere and the sublithospheric mantle "
			                       "folows dislocation creep and diffusion creep flow laws respectively (Karato and Wu, 1993)")
  }
}
