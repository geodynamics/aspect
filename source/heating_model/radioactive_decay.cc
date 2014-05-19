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
/*  $Id$  */


#include <aspect/heating_model/radioactive_decay.h>
#include <aspect/global.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    Radioactive_decay<dim>::Radioactive_decay ()
    {}



    template <int dim>
    double
    Radioactive_decay<dim>::
    specific_heating_rate (const double,
        const double,
        const std::vector<double> &,
        const Point<dim> &p) const
    {
        double timedependent_radioactive_heating_rate=0;
        if(num_radio_heating_elements!=0)
            for(unsigned i_radio=0;i_radio<num_radio_heating_elements;i_radio++)
                timedependent_radioactive_heating_rate+=
                    radioactive_heating_rate[i_radio]
                    *radioactive_initial_consentration[i_radio]*1e-6
                    *std::pow(0.5,time/half_decay_time[i_radio]);
        return (timedependent_radioactive_heating_rate);
    }


    template <int dim>
    void
    Radioactive_decay<dim>::update ()
    {
      time = this->get_time();
      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      if (this->convert_output_to_years())
        time/=year_in_seconds;
    }


    template <int dim>
    void
    Radioactive_decay<dim>::declare_parameters (ParameterHandler &prm)
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
                                "Half decay time. Units: (Seconds), or (Years) if set 'use years instead of seconds'.");
            prm.declare_entry("Initial consentration","",
                                Patterns::List (Patterns::Double (0)),
                                "Initial consentration of different elments (ppm)");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Radioactive_decay<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Radioactive decay");
        {
        
            num_radio_heating_elements= prm.get_integer ("Number of elements");
            radioactive_heating_rate=Utilities::string_to_double
                (Utilities::split_string_list
                (prm.get("Heating rate")));
            AssertThrow(radioactive_heating_rate.size()==num_radio_heating_elements,
                ExcMessage("Number of heating rate entities does not match "
                           "the number of radioactive elements."));
            radioactive_initial_consentration=Utilities::string_to_double
                (Utilities::split_string_list
                (prm.get("Initial consentration")));
            AssertThrow(radioactive_initial_consentration.size()==num_radio_heating_elements,
                ExcMessage("Number of initial consentration entities does not match "
                           "the number of radioactive elements."));
            half_decay_time=Utilities::string_to_double
                (Utilities::split_string_list
                 (prm.get("Half decay time")));
            AssertThrow(half_decay_time.size()==num_radio_heating_elements,
                ExcMessage("Number of half decay time entities does not match "
                           "the number of radioactive elements."));
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
    ASPECT_REGISTER_HEATING_MODEL(Radioactive_decay,
                                                 "radioactive decay",
                                                 "Implementation of a model in which the heating "
                                                 "rate is given in terms of an explicit formula "
                                                 "that is elaborated in the parameters in section "
                                                 "``Heating model|Function''. "
                                                 "\n\n"
                                                 "The formula is interpreted as having units "
                                                 "W/kg."
                                                 "\n\n"
                                                 "Since the symbol $t$ indicating time "
                                                 "may appear in the formulas for the heating rate"
                                                 ", it is interpreted as having units "
                                                 "seconds unless the global parameter ``Use "
                                                 "years in output instead of seconds'' is set.")
  }
}
