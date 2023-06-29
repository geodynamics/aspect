/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/termination_criteria/user_request.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    UserRequest<dim>::execute()
    {
      return Utilities::fexists(this->get_output_directory()+filename_to_test, this->get_mpi_communicator());
    }

    template <int dim>
    void
    UserRequest<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.enter_subsection("User request");
        {
          prm.declare_entry ("File name", "terminate-aspect",
                             Patterns::FileName(Patterns::FileName::input),
                             "The name of a file that, if it exists in the output "
                             "directory (whose name is also specified in the input file) "
                             "will lead to termination of the simulation. "
                             "The file's location is chosen to be in the output directory, "
                             "rather than in a generic location such as the ASPECT directory, "
                             "so that one can run multiple simulations at the same time (which "
                             "presumably write to different output directories) and can "
                             "selectively terminate a particular one.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    UserRequest<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Termination criteria");
      {
        prm.enter_subsection("User request");
        {
          filename_to_test = prm.get ("File name");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace TerminationCriteria
  {
    ASPECT_REGISTER_TERMINATION_CRITERION(UserRequest,
                                          "user request",
                                          "Terminate the simulation gracefully when a file with a specified "
                                          "name appears in the output directory. This allows the user to "
                                          "gracefully exit the simulation at any time by simply creating "
                                          "such a file using, for example, \\texttt{touch output/terminate}. "
                                          "The file's location is chosen to be in the output directory, "
                                          "rather than in a generic location such as the ASPECT directory, "
                                          "so that one can run multiple simulations at the same time (which "
                                          "presumably write to different output directories) and can "
                                          "selectively terminate a particular one.")
  }
}
