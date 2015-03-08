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


#include <aspect/postprocess/geoid.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace internal
    {
      HarmonicCoefficients::HarmonicCoefficients(const unsigned int max_degree)
      {
        //TODO: find real size
        sine_coefficients.resize(max_degree);
        cosine_coefficients.resize(max_degree);
      }

      template <int dim>
      SphericalHarmonicsExpansion<dim>::SphericalHarmonicsExpansion(const unsigned int max_degree)
      :
      coefficients(max_degree)
      {

        // Set up coefficients with right sizes according to max_degree
      }

      template <int dim>
      void
      SphericalHarmonicsExpansion<dim>::add_point (const Point<dim> &position,
                                                  const double value)
      {
        const std_cxx1x::array<dim,double> spherical_position =
            Utilities::spherical_coordinates(position);
        // calculate contribution of value to coefficients
      }

      template <int dim>
      HarmonicCoefficients
      SphericalHarmonicsExpansion<dim>::get_coefficients () const
      {
        // TODO: How to calculate the summed coefficients?
        // What about parallelization?

        return coefficients;
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    Geoid<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      const GeometryModel::SphericalShell<dim> *geometry_model = dynamic_cast<const GeometryModel::SphericalShell<dim> *>
                                                    (&this->get_geometry_model());

      AssertThrow (geometry_model != 0,
                   ExcMessage("The geoid postprocessor is currently only implemented for "
                       "the spherical shell geometry model."));

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

      // loop over all of the surface cells and if one less than h/3 away from
      // the top surface, evaluate the stress at its center
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
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
              // density and add its contribution to the spherical harmonics

              for (unsigned int q=0; q<quadrature_formula.size(); ++q)
                {
                  const Point<dim> location = fe_values.quadrature_point(q);
                  const double viscosity = out.viscosities[q];
                  const double density   = out.densities[q];
                  const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(location);

                  const Tensor<1,dim> gravity_direction = gravity/gravity.norm();

                  const unsigned int layer_id = static_cast<unsigned int> (geometry_model->depth(location) / geometry_model->maximal_depth() * (number_of_layers-1));

                  expansions[layer_id].AddPoint(location,density);
                }

              // see if the cell is at the *top* boundary
              bool surface_cell = false;
              for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
                if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
                  {
                    surface_cell = true;
                    break;
                  }

              if (surface_cell)
                {
                  // Add topography contribution
                }
            }

      //TODO: loop over all layer_ids and compute coefficients, afterwards sum coefficients for surface geoid
      HarmonicCoefficients geoid_expansion;
      for (unsigned int layer_id = 0; layer_id < number_of_layers; ++layer_id)
        {
          const HarmonicCoefficients layer_contribution = expansions[layer_id].get_coefficients();

          //TODO scaling for depth
          geoid_expansion += layer_contribution;
        }

      // Write the solution to an output file
      for (unsigned int i=0; i<stored_values.size(); ++i)
        {
          output << stored_values[i].first
                 << ' '
                 << stored_values[i].second -
                 (subtract_mean_dyn_topography
                  ?
                  average_topography
                  :
                  0.)
                 << std::endl;
        }


      const std::string filename = this->get_output_directory() +
                                   "dynamic_topography." +
                                   Utilities::int_to_string(this->get_timestep_number(), 5);

      const unsigned int max_data_length = Utilities::MPI::max (output.str().size()+1,
                                                                this->get_mpi_communicator());
      const unsigned int mpi_tag = 124;

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

    template <int dim>
    void
    Geoid<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          prm.declare_entry ("Include topography contribution", "false",
                             Patterns::Bool (),
                             "Option to include the contribution of dynamic "
                             "topography to the geoid.");
          prm.declare_entry ("Number of layers", "20",
                             Patterns::Integer (1),
                             "The geoid contribution is added on a per-layer basis. This parameter "
                             "sets the number of layers. Similar to the depth-average "
                             "postprocessor, the number of layers should correspond roughly to "
                             "the available model resolution.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Geoid<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Geoid");
        {
          include_topography_contribution   = prm.get_bool("Include topography contribution");
          number_of_layers                  = prm.get_integer("Number of layers");
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
    ASPECT_REGISTER_POSTPROCESSOR(Geoid,
                                  "geoid",
                                  "A postprocessor that computes a measure of geoid height "
                                  "based on the")
  }
}
