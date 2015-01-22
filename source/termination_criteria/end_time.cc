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


#include <aspect/termination_criteria/end_time.h>


namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    EndTime<dim>::execute()
    {
      return (this->get_time() > end_time);
    }

    template <int dim>
    double EndTime<dim>::check_for_last_time_step (const double time_step) const
    {
      if ((this->get_time()<end_time)
          &&
          (this->get_time()+time_step > end_time))
        return end_time - this->get_time();
      else
        return time_step;
    }

    template <int dim>
    void
    EndTime<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry ("End time",
                         /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                              year_in_seconds) = */ "5.69e+300",
                         Patterns::Double (),
                         "The end time of the simulation. The default value is a number "
                         "so that when converted from years to seconds it is approximately "
                         "equal to the largest number representable in floating point "
                         "arithmetic. For all practical purposes, this equals infinity. "
                         "Units: Years if the "
                         "'Use years in output instead of seconds' parameter is set; "
                         "seconds otherwise.");
    }


    template <int dim>
    void
    EndTime<dim>::parse_parameters (ParameterHandler &prm)
    {
      // read end time from parameter file. if it is to be interpreted
      // in years rather than seconds, then do the conversion
      end_time = prm.get_double ("End time");
      if (prm.get_bool ("Use years in output instead of seconds") == true)
        end_time *= year_in_seconds;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TerminationCriteria
  {
    ASPECT_REGISTER_TERMINATION_CRITERION(EndTime,
                                          "end time",
                                          "Terminate the simulation once the end time "
                                          "specified in the input file has been reached. "
                                          "Unlike all other termination criteria, this "
                                          "criterion is \\textit{always} active, whether it "
                                          "has been explicitly selected or not in the input file "
                                          "(this is done to preserve historical behavior of "
                                          "\\aspect{}, but it also likely does not inconvenience "
                                          "anyone since it is what would be selected in most "
                                          "cases anyway).")
  }
}
