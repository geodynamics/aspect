// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
          AssertThrow(std::system(nullptr) != 0,
                      ExcMessage("The \"command\" postprocessor required a command-processor, "
                                 "which appears to be unavailable on this system."));

          const int error = std::system(command.c_str());

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
