/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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

#include <aspect/material_model/perplex_simpler.h>

extern "C" { 
#include <perplex_c.h>
}

namespace aspect
{
  namespace MaterialModel
  {
    
    template <int dim>
    void
    PerpleXSimpler<dim>::initialize()
    {
      ini_phaseq(perplex_file_name.c_str()); // this line initializes meemum
      
      // Initialize the various arrays
      wtphases = new double[p_size_phases];
      cphases = new double[p_size_phases * p_size_components];
      sysprop = new double[p_size_sysprops]; 
      namephases = new char[p_size_phases * p_pname_len];
	
    }
    
    template <int dim>
    bool
    PerpleXSimpler<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    double
    PerpleXSimpler<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    void
    PerpleXSimpler<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      /* Instead of evaluating at every quadrature point per cell,
       * we here average the P, T and X values, and evaluate once.
       * This is much quicker than evaluating at all quadrature 
       * points, and if the grid is fine, it should be a reasonable
       * approximation
       */

      unsigned int n_quad = in.position.size(); // number of quadrature points in cell
      unsigned int n_comp = in.composition[0].size(); // number of components in rock
      
      const double average_temperature = std::min(max_temperature,
						  std::max(min_temperature,
							   (accumulate( in.temperature.begin(), in.temperature.end(), 0.0) /
							    n_quad)));
      const double average_pressure = std::min(max_pressure,
						  std::max(min_pressure,
							   (accumulate( in.pressure.begin(), in.pressure.end(), 0.0) /
							    n_quad)));
      
      std::vector<double> comp;
      comp.resize(n_comp);
      
      for (unsigned int c=0; c<n_comp; ++c)
	{
	  for (unsigned int i=0; i<n_quad; ++i)
	    {   
	      comp[c] += in.composition[i][c];
	      out.reaction_terms[i][c] = 0.0;
	    }
	  comp[c] /= (double)n_quad;
	}

      // Here is the call to PerpleX/meemum
      int nphases;
      phaseq(average_pressure/1.e5, average_temperature,
	     n_comp, comp.data(), &nphases, wtphases, cphases,
	     sysprop, namephases, phaseq_dbg);

      const std::string condition = std::to_string(average_pressure) +" bar, " + std::to_string(average_temperature) + " K";
      AssertThrow(isnan(sysprop[9]) == false, ExcMessage("PerpleX returned NaN for density at " + condition + ". Aborting. "
							 "Please adjust the P-T bounds in the parameter file or adjust the PerpleX files."));
      AssertThrow(isnan(sysprop[11]) == false, ExcMessage("PerpleX returned NaN for heat capacity at " + condition + ". Aborting. "
							  "Please adjust the P-T bounds in the parameter file or adjust the PerpleX files."));
      AssertThrow(isnan(sysprop[12]) == false, ExcMessage("PerpleX returned NaN for thermal expansivity at " + condition + ". Aborting. "
							  "Please adjust the P-T bounds in the parameter file or adjust the PerpleX files."));
      AssertThrow(isnan(sysprop[13]) == false, ExcMessage("PerpleX returned NaN for compressibility at " + condition + ". Aborting. "
							  "Please adjust the P-T bounds in the parameter file or adjust the PerpleX files."));
      
      
      for (unsigned int i=0; i<n_quad; ++i)
        {   
          out.viscosities[i] = eta;
          out.thermal_conductivities[i] = k_value;
	  out.densities[i] = sysprop[9];
	  out.specific_heat[i] = sysprop[11]*(1000./sysprop[16]); // molar Cp * (1000/molar mass) (g)
	  out.thermal_expansion_coefficients[i] = sysprop[12];
	  out.compressibilities[i] = sysprop[13]*1.e5;
        }

    }


    template <int dim>
    void
    PerpleXSimpler<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("PerpleX simpler model");
        {
	  
          prm.declare_entry ("PerpleX input file name", "rock.dat",
                             Patterns::Anything (),
                             "The name of the PerpleX input file (should end with .dat).");
          prm.declare_entry ("Viscosity", "5e24",
                             Patterns::Double (0),
                             "The value of the viscosity $\\eta$. Units: $kg/m/s$.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: $W/m/K$.");
          prm.declare_entry ("Minimum material temperature", "0.",
                             Patterns::Double (0),
                             "The value of the minimum temperature used to query PerpleX. "
                             "Units: $K$.");
          prm.declare_entry ("Maximum material temperature", "6000.",
                             Patterns::Double (0),
                             "The value of the maximum temperature used to query PerpleX. "
                             "Units: $K$.");
          prm.declare_entry ("Minimum material pressure", "1.e5",
                             Patterns::Double (0),
                             "The value of the minimum pressure used to query PerpleX. "
                             "Units: $Pa$.");
          prm.declare_entry ("Maximum material pressure", "1.e12",
                             Patterns::Double (0),
                             "The value of the maximum pressure used to query PerpleX. "
                             "Units: $Pa$.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    PerpleXSimpler<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("PerpleX simpler model");
        {
	  perplex_file_name   = prm.get ("PerpleX input file name");
          eta                 = prm.get_double ("Viscosity");
          k_value             = prm.get_double ("Thermal conductivity");
          min_temperature     = prm.get_double ("Minimum material temperature");
          max_temperature     = prm.get_double ("Maximum material temperature");
          min_pressure        = prm.get_double ("Minimum material pressure");
          max_pressure        = prm.get_double ("Maximum material pressure");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      
      this->model_dependence.density = NonlinearDependence::temperature
	                               | NonlinearDependence::pressure
	                               | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::temperature
                                               | NonlinearDependence::pressure
                                               | NonlinearDependence::compositional_fields;
      this->model_dependence.specific_heat = NonlinearDependence::temperature
                                             | NonlinearDependence::pressure
                                             | NonlinearDependence::compositional_fields;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PerpleXSimpler,
                                   "perplex simpler",
                                   "A material model that has constant values "
                                   "for viscosity and thermal conductivity, and "
				   "calculates other properties on-the-fly using "
				   "PerpleX meemum. Compositional fields correspond "
				   "to the individual components in the order given "
				   "in the PerpleX file.")
  }
}
