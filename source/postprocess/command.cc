/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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

#include <aspect/postprocess/command.h>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    Command<dim>::execute (TableHandler &)
    {
      if (on_all_processes ||
          (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0))
        {

          // Check if a command-processor is available by calling system() with a
          // null pointer. System is guaranteed to return non-zero if it finds
          // a terminal and zero if there is none (like on the compute nodes of
          // some cluster architectures, e.g. IBM BlueGene/Q)
          AssertThrow(system((char *)0) != 0,
                      ExcMessage("The \"command\" postprocessor required a command-processor, "
                                 "which appears to be unavailable on this system."));

          const int error = system(command.c_str());

          if (error != 0)
            {
              std::string err_str = (error<0 ? "-" : "") +
                                    Utilities::int_to_string(std::abs(error));


              if (terminate_on_failure)
                {
                  AssertThrow(false, ExcMessage("Command <" +
                                                command +
                                                "> failed with error code: " +
                                                err_str +
                                                ". Terminating process."));
                }
              else
                {
                  std::cerr << "*** WARNING: Command <" << command
                            << "> failed with error code: "
                            << err_str
                            << ". Continuing anyway."
                            << std::endl;
                }
            }
        }

      return std::make_pair (std::string("Running command:"),
                             command);
    }

    template <int dim>
    void
    Command<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Command");
        {
          prm.declare_entry ("Terminate on failure", "false",
                             Patterns::Bool(),
                             "Select whether \\aspect{} should terminate if the "
                             "command returns a non-zero exit status.");

          prm.declare_entry ("Run on all processes", "false",
                             Patterns::Bool(),
                             "Whether to run command from all processes (true), "
                             "or only on process 0 (false).");

          prm.declare_entry ("Command", "",
                             Patterns::Anything(),
                             "Command to execute.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Command<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Command");
        {
          terminate_on_failure = prm.get_bool ("Terminate on failure");
          on_all_processes = prm.get_bool ("Run on all processes");
          command = prm.get ("Command");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Command,
                                  "command",
                                  "A postprocessor that executes a command line "
                                  "process.")
  }
}
