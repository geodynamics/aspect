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
/*  $Id$  */


#include <aspect/termination_criteria/end_time.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    EndTime<dim>::execute()
    {
      return (this->get_time() >= end_time);
    }

    template <int dim>
    void
    EndTime<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry ("End time", "1e300",
                         Patterns::Double (),
                         "The end time of the simulation. Units: years if the "
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
      if (this->convert_output_to_years() == true)
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
                                          "end_time",
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
