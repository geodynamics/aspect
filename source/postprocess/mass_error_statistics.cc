/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/mass_error_statistics.h>
#include <aspect/material_model/interface.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MassErrorStatistics<dim>::execute (TableHandler &statistics)
    {
      const FE_Q<dim> density_element(this->get_fe().base_element(this->introspection().base_elements.velocities).degree);
      const std::vector<Point<dim> > support_points = density_element.get_unit_support_points();

      const QGauss<dim> quadrature_formula (density_element.degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      /* We need three different FEValues objects for the operations below.
       * fe_values_support evaluates the material model at the support points
       * of the hypothetical density element declared above.
       * fe_values directly evaluates the velocity at the quadrature points
       * fe_values_density is used to interpolate the density from the support
       * points to the quadrature points.
       */

      FEValues<dim> fe_values_support (this->get_mapping(),
                                       this->get_fe(),
                                       Quadrature<dim>(support_points),
                                       update_q_points |
                                       update_values |
                                       update_gradients);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values |
                               update_gradients |
                               update_JxW_values |
                               update_q_points);

      FEValues<dim> fe_values_density (this->get_mapping(),
                                       density_element,
                                       quadrature_formula,
                                       update_values |
                                       update_gradients);


      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      MaterialModel::MaterialModelInputs<dim> in(support_points.size(), this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(support_points.size(), this->n_compositional_fields());

      double local_integral_mass_error = 0;
      double local_norm_mass_error = 0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values_support.reinit (cell);
            fe_values_density.reinit(typename Triangulation<dim>::active_cell_iterator(cell));

            fe_values_support[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                in.velocity);
            fe_values_support[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                in.temperature);
            fe_values_support[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                in.pressure);
            fe_values_support[this->introspection().extractors.pressure].get_function_gradients (this->get_solution(),
                in.pressure_gradient);
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values_support[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                  composition_values[c]);

            in.position = fe_values_support.get_quadrature_points();

            // since we are not reading the viscosity and the viscosity
            // is the only coefficient that depends on the strain rate,
            // we need not compute the strain rate. set the corresponding
            // array to empty, to prevent accidental use and skip the
            // evaluation of the strain rate in evaluate().
            in.strain_rate.resize(0);

            for (unsigned int i=0; i<support_points.size(); ++i)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[i][c] = composition_values[c][i];
              }
            in.cell = &cell;

            this->get_material_model().evaluate(in, out);

            std::vector<Tensor<1,dim> > velocity(n_q_points);
            std::vector<double>         velocity_divergence(n_q_points);
            fe_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(),velocity);
            fe_values[this->introspection().extractors.velocities].get_function_divergences(this->get_solution(),velocity_divergence);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                Tensor<1,dim> density_gradient = Tensor<1,dim>();
                double        density = 0.0;
                for (unsigned int i = 0; i < support_points.size(); ++i)
                  {
                    density_gradient += out.densities[i] * fe_values_density.shape_grad(i,q);
                    density += out.densities[i] * fe_values_density.shape_value(i,q);
                  }

                const double point_error = (density * velocity_divergence[q]
                                            + velocity[q] * density_gradient);

                local_integral_mass_error += point_error * fe_values.JxW(q);
                local_norm_mass_error += point_error * point_error * fe_values.JxW(q);
              }
          }

      const double global_integral_mass_error
        = Utilities::MPI::sum (local_integral_mass_error, this->get_mpi_communicator());
      const double global_norm_mass_error
        = std::sqrt(Utilities::MPI::sum (local_norm_mass_error, this->get_mpi_communicator()));

      if (this->convert_output_to_years() == true)
        {
          statistics.add_value ("Integral mass error (kg/year)",
                                global_integral_mass_error * year_in_seconds);
          statistics.add_value ("Norm mass error (kg/year)",
                                global_norm_mass_error * year_in_seconds);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "Integral mass error (kg/year)",
                                      "Norm mass error (kg/year)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }
      else
        {
          statistics.add_value ("Integral mass error (kg/s)",
                                global_integral_mass_error);
          statistics.add_value ("Norm mass error (kg/s)",
                                global_norm_mass_error);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "Integral mass error (kg/s)",
                                      "Norm mass error (kg/s)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << global_integral_mass_error *year_in_seconds
               << " kg/year, "
               << global_norm_mass_error *year_in_seconds
               << " kg/year";
      else
        output << global_integral_mass_error
               << " kg/s, "
               << global_norm_mass_error
               << " kg/s";

      return std::pair<std::string, std::string> ("Integral mass error, Norm mass error:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MassErrorStatistics,
                                  "mass error statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "error of the mass conservation equation. More precisely, it "
                                  "computes the norm and the integral of the equation "
                                  "$\\rho \\Nabla \\mathbf u + \\Nabla \\rho \\cdot \\mathbf u$. "
                                  "Therefore, its results indicate how accurate the chosen approximation "
                                  "for compressibility is for the given model setup.")
  }
}
