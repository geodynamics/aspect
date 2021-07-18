/*
  Copyright (C) 2021 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/postprocess/adjoint_kernels.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    std::pair<std::string,std::string>
    AdjointKernels<dim>::execute (TableHandler &)
    {
      if (this->n_compositional_fields() == 0)
        return std::pair<std::string,std::string>();

      // This postprocessor is only useful when solving the adjoint Stokes equations so throw
      // an exception if a different solver is chosen.
      Assert(this->get_parameters().nonlinear_solver == Parameters<dim>::NonlinearSolver::no_Advection_adjoint_Stokes,
             ExcMessage("Error: Postprocessor Adjoint kernels is only sensible with the no Advection, adjoint Stokes solver scheme."))

      // If the no Advection, adjoint Stokes solver scheme is chosen the parameter file is already checked such
      // that it has at least 2 compositional fields and that the compositional fields have the correct names.

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);

      const unsigned int n_q_points = quadrature_formula.size();
      std::vector<double> compositional_values_density(n_q_points);
      std::vector<double> compositional_values_viscosity(n_q_points);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      std::ostringstream output;
      std::vector<double> kernel_density_vec;
      std::vector<double> kernel_viscosity_vec;
      std::vector<Point<dim> > location;

      // Prevent frequent memory reallocation when the vector grows
      kernel_density_vec.reserve(this->get_triangulation().n_locally_owned_active_cells());
      kernel_viscosity_vec.reserve(this->get_triangulation().n_locally_owned_active_cells());
      location.reserve(this->get_triangulation().n_locally_owned_active_cells());

      // get the index for the two compositional fields that correspond to the kernels
      const unsigned int density_idx = this->introspection().compositional_index_for_name("density_increment");
      const unsigned int viscosity_idx = this->introspection().compositional_index_for_name("viscosity_increment");

      // compute the integral quantities by quadrature
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            double kernel_density = 0;
            double kernel_viscosity = 0;
            double volume = 0;
            const bool use_manifold_in_curved_meshes = true;
            const Point<dim> midpoint_of_cell = cell->center(use_manifold_in_curved_meshes);

            fe_values[this->introspection().extractors.compositional_fields[density_idx]].get_function_values (this->get_solution(),
                compositional_values_density);
            fe_values[this->introspection().extractors.compositional_fields[viscosity_idx]].get_function_values (this->get_solution(),
                compositional_values_viscosity);

            for (unsigned int q=0; q<n_q_points; ++q)
              {
                kernel_density += compositional_values_density[q]*fe_values.JxW(q);
                kernel_viscosity += compositional_values_viscosity[q]*fe_values.JxW(q);
                volume += fe_values.JxW(q);
              }

            // safe cell value into vector for print
            kernel_density_vec.push_back(kernel_density/volume);
            kernel_viscosity_vec.push_back(kernel_viscosity/volume);
            location.push_back(midpoint_of_cell);
          }

      // Write the solution to an output file
      for (unsigned int i=0; i<kernel_viscosity_vec.size(); ++i)
        {
          output << location[i]
                 << ' '
                 << kernel_density_vec[i]
                 << ' '
                 << kernel_viscosity_vec[i]
                 << std::endl;
        }

      std::string filename = this->get_output_directory() +
                             "adjoint_kernel." +
                             Utilities::int_to_string(this->get_timestep_number(), 5);
      if (this->get_parameters().run_postprocessors_on_nonlinear_iterations)
        filename.append("." + Utilities::int_to_string (this->get_nonlinear_iteration(), 4));

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 123;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

          file << "# "
               << ((dim==2)? "x y" : "x y z")
               << " density  viscosity " << std::endl;

          // first write out the data we have created locally
          file << output.str();

          std::string tmp;
          tmp.resize (max_data_length, '\0');

          // then loop through all of the other processors and collect
          // data, then write it to the file
          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              // get the data. note that MPI says that an MPI_Recv may receive
              // less data than the length specified here. since we have already
              // determined the maximal message length, we use this feature here
              // rather than trying to find out the exact message length with
              // a call to MPI_Probe.
              MPI_Recv (&tmp[0], max_data_length, MPI_CHAR, p, mpi_tag,
                        this->get_mpi_communicator(), &status);

              // output the string. note that 'tmp' has length max_data_length,
              // but we only wrote a certain piece of it in the MPI_Recv, ended
              // by a \0 character. write only this part by outputting it as a
              // C string object, rather than as a std::string
              file << tmp.c_str();
            }
        }
      else
        // on other processors, send the data to processor zero. include the \0
        // character at the end of the string
        {
          MPI_Send (&output.str()[0], output.str().size()+1, MPI_CHAR, 0, mpi_tag,
                    this->get_mpi_communicator());
        }

      return std::pair<std::string,std::string> ("Writing adjoint kernels:", filename);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(AdjointKernels,
                                  "adjoint kernels",
                                  "A postprocessor that writes out the sensitivity kernel for "
                                  "density and viscosity. These sensitivity kernels are calculated "
                                  "based on the solution of the forward and adjoint Stokes equations. "
                                  "This postprocessor is therefore only sensible if the solver "
                                  "scheme is no Advection, adjoint Stokes. The sensitivity "
                                  "kernels are computed by the helper function compute_sensitivity_kernel "
                                  "and are stored in the two compositional fields. Consequently this "
                                  "postprocessor requires at least 2 compositional fields (to store "
                                  "the two kernels) and the fields have to have the names density_increment "
                                  "and viscosity_increment."
                                  "\n"
                                  "The file format consists of lines with Euclidean coordinates "
                                  "followed by the corresponding density and then viscosity sensitivity "
                                  "kernel value."
                                  "\n\n")
  }
}
