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


#include <aspect/termination_criteria/user_request.h>

namespace aspect
{
  namespace TerminationCriteria
  {
    template <int dim>
    bool
    UserRequest<dim>::execute(void)
    {
      std::fstream        check_file;

      // Only check for the file on the root process to avoid overloading the filesystem
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          check_file.open((this->get_output_directory()+_end_filename).c_str());
          if (check_file.is_open())
            {
              check_file.close();
              return true;
            }
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
          prm.declare_entry ("File name", "terminate_aspect",
                             Patterns::FileName(Patterns::FileName::input),
                             "");
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
          _end_filename = prm.get ("File name");
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
                                          "user_request",
                                          "Terminate the simulation gracefully when a file with a specified "
                                          "name appears in the output directory. This allows the user to "
                                          "gracefully exit the simulation at any time.")
  }
}
