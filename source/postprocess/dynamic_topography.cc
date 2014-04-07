/*
  Copyright (C) 2011, 2012, 2013, 2014 by the authors of the ASPECT code.

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


#include <aspect/postprocess/dynamic_topography.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DynamicTopography<dim>::execute (TableHandler &statistics)
    {
      const QMidpoint<dim> quadrature_formula;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_q_points);

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(fe_values.n_quadrature_points, this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(fe_values.n_quadrature_points, this->n_compositional_fields());

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      // have a stream into which we write the data. the text stream is then
      // later sent to processor 0
      std::ostringstream output;

      // loop over all of the surface cells and if one less than h/2 away from
      // the top surface, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          if (cell->at_boundary() &&
              (this->get_geometry_model().depth (cell->center()) < cell->diameter()/2))
            {
              fe_values.reinit (cell);

              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature]
              .get_function_values (this->get_solution(), in.temperature);
              fe_values[this->introspection().extractors.pressure]
              .get_function_values (this->get_solution(), in.pressure);
              fe_values[this->introspection().extractors.velocities]
              .get_function_symmetric_gradients (this->get_solution(), in.strain_rate);

              in.position = fe_values.get_quadrature_points();

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]]
                .get_function_values(this->get_solution(),
                                     composition_values[c]);
              for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                {
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    in.composition[i][c] = composition_values[c][i];
                }

              this->get_material_model().evaluate(in, out);


              // for each of the quadrature points, evaluate the
              // stress and compute the component in direction of the
              // gravity vector
              for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                {
                  const Point<dim> location = fe_values.quadrature_point(q);
                  const double viscosity = out.viscosities[q];
                  const double density   = out.densities[q];

//TODO: We need to subtract 2/3*div(u) from the stress here in the compressible case
                  const SymmetricTensor<2,dim> stress = 2 * viscosity * in.strain_rate[q];

                  const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);
                  const Tensor<1,dim> gravity_direction = gravity/gravity.norm();

                  const double sigma_rr           = gravity_direction * (stress * gravity_direction);
                  const double dynamic_topography = sigma_rr / gravity.norm() / density;

                  output << location
                         << ' '
                         << dynamic_topography
                         << std::endl;
                }
            }

      const std::string filename = this->get_output_directory() +
                                   "dynamic_topography." +
                                   Utilities::int_to_string(this->get_timestep_number(), 5);

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 123;

      // on processor 0, collect all of the data the individual processors send
      // and concatenate them into one file
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::ofstream file (filename.c_str());

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

      return std::pair<std::string,std::string>("Writing dynamic topography:",
                                                filename);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DynamicTopography,
                                  "dynamic topography",
                                  "A postprocessor that computes a measure of dynamic topography "
                                  "based on the stress at the surface. The data is writte into a "
                                  "file named 'dynamic_topography.NNNNN' in the output directory, "
                                  "where NNNNN is the number of the time step."
                                  "\n\n"
                                  "The exact approach works as follows: At the centers of all cells "
                                  "that sit along the top surface, we evaluate the stress and "
                                  "evaluate the component of it in the direction in which "
                                  "gravity acts. In other words, we compute "
                                  "$\\sigma_{rr}={\\hat g}^T(2 \\eta \\varepsilon(\\mathbf u))\\hat g$ "
                                  "where $\\hat g = \\mathbf g/\\|\\mathbf g\\|$ is the direction of "
                                  "the gravity vector $\\mathbf g$. From this, the dynamic "
                                  "topography is computed using the formula "
                                  "$h=\\frac{\\sigma_{rr}}{\\|\\mathbf g\\| \\rho}$ where $\\rho$ "
                                  "is the density at the cell center."
                                  "\n\n"
                                  "(As a side note, the postprocessor chooses the cell center "
                                  "instead of the center of the cell face at the surface, where we "
                                  "really are interested in the quantity, since "
                                  "this often gives better accuracy. The results should in essence "
                                  "be the same, though.)")
  }
}
