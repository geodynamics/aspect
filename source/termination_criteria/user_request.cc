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


#include <aspect/termination_criteria/user_request.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    UserRequest<dim>::execute(void)
    {
      // Only check for the file on the root process to avoid overloading the filesystem.
      // The plugin manager later does an OR operation over all
      // processors
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::fstream check_file((this->get_output_directory()+filename_to_test).c_str());
          return check_file.is_open();
        }
      return false;
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
                             "rather than in a generic location such as the Aspect directory, "
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
                                          "rather than in a generic location such as the Aspect directory, "
                                          "so that one can run multiple simulations at the same time (which "
                                          "presumably write to different output directories) and can "
                                          "selectively terminate a particular one.")
  }
}
