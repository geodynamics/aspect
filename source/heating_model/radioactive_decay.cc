/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#include <aspect/heating_model/radioactive_decay.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/global.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    RadioactiveDecay<dim>::RadioactiveDecay ()
    {}



    template <int dim>
    double
    RadioactiveDecay<dim>::
    specific_heating_rate (const double,
        const double,
        const std::vector<double> &composition,
        const Point<dim> &position) const
    {
        double timedependent_radioactive_heating_rate=0;
        if(n_radio_heating_elements!=0)
        {
            double crust_percent=0;
            if(is_crust_defined_by_composition)
            {
                AssertThrow(crust_composition_num < composition.size(), ExcMessage("The composition number of crust is "
                                " larger than number of composition fields."));
                crust_percent=composition[crust_composition_num];
                if(crust_percent<0)crust_percent=0;
                if(crust_percent>1)crust_percent=1;
            }
            else
                if((this->get_geometry_model()).depth(position) < crust_depth)
                    crust_percent=1.;
                    
            for(unsigned i_radio=0;i_radio<n_radio_heating_elements;i_radio++)
                timedependent_radioactive_heating_rate+=
                    radioactive_heating_rate[i_radio]
                    *(radioactive_initial_concentration_mantle[i_radio]*(1-crust_percent)
                    +radioactive_initial_concentration_crust[i_radio]*crust_percent)*1e-6
                    *std::pow(0.5,time/half_decay_time[i_radio]);
        }
        return (timedependent_radioactive_heating_rate);
    }


    template <int dim>
    void
    RadioactiveDecay<dim>::update ()
    {
      time = this->get_time();
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        time/=year_in_seconds;
    }


    template <int dim>
    void
    RadioactiveDecay<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Radioactive decay");
        {
            prm.declare_entry("Number of elements","0",
                                Patterns::Integer(0),
                                "Number of radioactive elements");
            prm.declare_entry("Heating rate","",
                                Patterns::List (Patterns::Double ()),
                                "Heating rate of different element (W/kg)");
            prm.declare_entry("Half decay time","",
                                Patterns::List (Patterns::Double (0)),
                                "Half decay time. Units: (Seconds), or "
                                "(Years) if set 'use years instead of seconds'.");
            prm.declare_entry("Initial concentration crust","",
                                Patterns::List (Patterns::Double (0)),
                                "Initial concentration of different elements (ppm)");
            prm.declare_entry("Initial concentration mantle","",
                                Patterns::List (Patterns::Double (0)),
                                "Initial concentration of different elements (ppm)");
            prm.declare_entry("Crust defined by composition","false",
                                Patterns::Bool(),
                                "Whether crust defined by composition or depth");
            prm.declare_entry("Crust depth","0",
                                Patterns::Double(),
                                "Depth of the crust when crust if defined by depth. "
                                "Units: m");
            prm.declare_entry("Crust composition number","0",
                                Patterns::Integer(0),
                                "Which composition field should be treated as crust");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    RadioactiveDecay<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Radioactive decay");
        {
        
            n_radio_heating_elements= prm.get_integer ("Number of elements");
            radioactive_heating_rate=Utilities::string_to_double
                (Utilities::split_string_list
                (prm.get("Heating rate")));
            AssertThrow(radioactive_heating_rate.size()==n_radio_heating_elements,
                ExcMessage("Number of heating rate entities does not match "
                           "the number of radioactive elements."));
                           
            half_decay_time=Utilities::string_to_double
                (Utilities::split_string_list
                 (prm.get("Half decay time")));
            AssertThrow(half_decay_time.size()==n_radio_heating_elements,
                ExcMessage("Number of half decay time entities does not match "
                           "the number of radioactive elements."));
                           
            radioactive_initial_concentration_crust=Utilities::string_to_double
                (Utilities::split_string_list
                (prm.get("Initial concentration crust")));
            AssertThrow(radioactive_initial_concentration_crust.size()==n_radio_heating_elements,
                        ExcMessage("Number of initial concentration entities in crust "
                                   "does not match the number of radioactive elements."));
                           
            radioactive_initial_concentration_mantle=Utilities::string_to_double
                (Utilities::split_string_list
                (prm.get("Initial concentration mantle")));
            AssertThrow(radioactive_initial_concentration_mantle.size()==n_radio_heating_elements,
                        ExcMessage("Number of initial concentration entities in mantle "
                                   "does not match the number of radioactive elements."));

            is_crust_defined_by_composition = prm.get_bool    ("Crust defined by composition");
            crust_depth                     = prm.get_double  ("Crust depth");
            crust_composition_num           = prm.get_integer ("Crust composition number");
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
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(RadioactiveDecay,
                                  "radioactive decay",
                                  "Implementation of a model in which the internal "
                                  "heating rate is decaying exponentially over time \n"
                                  "The formula is interpreted as having units W/kg. \n\n")
  }                                              
}                                                
                                                 
                                                                