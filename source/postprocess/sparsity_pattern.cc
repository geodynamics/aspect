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

#include <aspect/postprocess/sparsity_pattern.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <fstream>
#include <iostream>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    SparsityPattern<dim>::SparsityPattern ()
      :
      // The following value is later read from the input file
      output_interval (0),
      // Initialize this to a nonsensical value; set it to the actual time
      // the first time around we get to check it
      next_output_time (std::numeric_limits<double>::quiet_NaN())
    {}

    template <int dim>
    std::pair<std::string,std::string>
    SparsityPattern<dim>::execute (TableHandler &statistics)
    {
      // If this is the first time we get here, set the next output time
      // to the current time. This makes sure we always produce data during
      // the first time step
      if (std::isnan(next_output_time))
        next_output_time = this->get_time();

      static int file_index;

      // See if output is requested at this time
      if (this->get_time() < next_output_time)
        return std::pair<std::string,std::string>();

      // If we're here, then output is requested, and we need to generate a new
      // sparsity pattern output file
      std::string data_file_name_prefix = this->get_output_directory()  + "Sparsity-" + Utilities::int_to_string(file_index++, 5);
      std::string data_file_name = data_file_name_prefix + (Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) > 0 ? "-" + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()),5) : "");
      // TODO: make this use file_index like the particle\output. This should be unique enough
      //       for now, but is a terrible solution
      std::ofstream output (data_file_name.c_str());
      // error handling
      if (!output) std::cout << "ERROR: could not create " << data_file_name << std::endl;

      // Get the DoFHandler from the simulation in order do calculate the
      // sparsity pattern
      const DoFHandler<dim> &dof_handler = this->get_dof_handler();

      CompressedSparsityPattern compressed_sparsity_pattern(dof_handler.n_dofs(),
                                                            dof_handler.n_dofs());

      DoFTools::make_sparsity_pattern (dof_handler, compressed_sparsity_pattern);

      dealii::SparsityPattern sparsity_pattern;
      sparsity_pattern.copy_from (compressed_sparsity_pattern);

      // Output the sparsity pattern in a format to be read by gnuplot
      // TODO: don't use gnuplot
      sparsity_pattern.print_gnuplot (output);

      output.close();

      // Set the next output time
      set_next_output_time(this->get_time());

      // Return the filename that was created
      return std::make_pair (std::string ("Writing sparsity pattern"),
                             data_file_name_prefix);
    }

    template <int dim>
    void
    SparsityPattern<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("sparsity pattern");
        {
          prm.declare_entry ("Time between output", "1e8",
                             Patterns::Double (0),
                             "The time interval between each generation of "
                             "sparsity pattern files. A value of zero indicates "
                             "that output should be generated in each time step. "
                             "Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    SparsityPattern<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Sparsity pattern");
        {
          output_interval = prm.get_double ("Time between output");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    SparsityPattern<dim>::set_next_output_time (const double current_time)
    {
      // if output_interval is positive, then set the next output interval to
      // a positive multiple.
      if (output_interval > 0)
        {
          // the current time is always in seconds, so we need to convert the output_interval to the same unit
          double output_interval_in_s = (this->convert_output_to_years()) ? (output_interval*year_in_seconds) : output_interval;

          // we need to compute the smallest integer that is bigger than current_time/my_output_interval,
          // even if it is a whole number already (otherwise we output twice in a row)
          next_output_time = (std::floor(current_time/output_interval_in_s)+1.0) * output_interval_in_s;
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(SparsityPattern,
                                  "sparsity pattern",
                                  "A postprocessor that generates a plot of the sparsity pattern of the finite element matrix")
  }
}



