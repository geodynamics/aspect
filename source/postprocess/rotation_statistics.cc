/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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


#include <aspect/postprocess/rotation_statistics.h>
#include <aspect/material_model/simple.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      void add_scientific_column(const std::string &name,
                                 const double value,
                                 TableHandler &statistics)
      {
        statistics.add_value(name,value);
        statistics.set_precision (name, 8);
        statistics.set_scientific (name, true);
      }
    }

    template <int dim>
    std::pair<std::string,std::string>
    RotationStatistics<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities).degree+1);

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, this->n_compositional_fields());

      Tensor<1,dim> local_angular_momentum;
      double local_scalar_angular_momentum = 0.0;

      SymmetricTensor<2,dim> local_moment_of_inertia;
      double local_scalar_moment_of_inertia = 0.0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // Set use_strain_rates to false since we don't need viscosity
            in.reinit(fe_values, cell, this->introspection(), this->get_solution(), false);
            if (!use_constant_density)
              {
                this->get_material_model().evaluate(in, out);
              }

            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q)
              {
                const double rho = (use_constant_density)
                                   ?
                                   1.0
                                   :
                                   out.densities[q];

                if (dim == 2)
                  {
                    // Get the velocity perpendicular to the position vector
                    const Tensor<1,dim> r_perp = cross_product_2d(in.position[q]);
                    // calculate a signed scalar angular momentum
                    local_scalar_angular_momentum += in.velocity[q] * r_perp * rho * fe_values.JxW(q);
                    local_scalar_moment_of_inertia += in.position[q].norm_square() * rho * fe_values.JxW(q);
                  }
                else if (dim == 3)
                  {
                    const Point<dim> position = in.position[q];
                    const Tensor<1,dim> r_cross_v = cross_product_3d(in.position[q], in.velocity[q]);

                    local_angular_momentum += r_cross_v * rho * fe_values.JxW(q);

                    // calculate moment of inertia
                    local_moment_of_inertia[0][0] += (position[1] * position[1] + position[2] * position[2]) * rho * fe_values.JxW(q);
                    local_moment_of_inertia[1][1] += (position[0] * position[0] + position[2] * position[2]) * rho * fe_values.JxW(q);
                    local_moment_of_inertia[2][2] += (position[0] * position[0] + position[1] * position[1]) * rho * fe_values.JxW(q);
                    local_moment_of_inertia[0][1] -= (position[0] * position[1]) * rho * fe_values.JxW(q);
                    local_moment_of_inertia[0][2] -= (position[0] * position[2]) * rho * fe_values.JxW(q);
                    local_moment_of_inertia[1][2] -= (position[1] * position[2]) * rho * fe_values.JxW(q);
                  }
              }
          }

      double global_angular_velocity = 0.0;
      double global_angular_momentum = 0.0;
      double global_moment_of_inertia = 0.0;
      SymmetricTensor<2,dim> global_tensor_moment_of_inertia;

      if (dim == 2)
        {
          global_angular_momentum =
            Utilities::MPI::sum (local_scalar_angular_momentum, this->get_mpi_communicator());
          global_moment_of_inertia =
            Utilities::MPI::sum (local_scalar_moment_of_inertia, this->get_mpi_communicator());
          global_angular_velocity = global_angular_momentum / global_moment_of_inertia;
        }
      else
        {
          const Tensor<1,dim> global_tensor_angular_momentum =
            Utilities::MPI::sum (local_angular_momentum, this->get_mpi_communicator());
          global_tensor_moment_of_inertia =
            Utilities::MPI::sum (local_moment_of_inertia, this->get_mpi_communicator());

          const SymmetricTensor<2,dim> inverse_moment(invert(Tensor<2,dim>(global_tensor_moment_of_inertia)));
          const Tensor<1,dim> global_vector_angular_velocity = inverse_moment * global_tensor_angular_momentum;

          global_angular_velocity = global_vector_angular_velocity.norm();
          global_angular_momentum = global_tensor_angular_momentum.norm();

          // Compute the moment of inertia around the current rotation axis
          global_moment_of_inertia = (global_tensor_moment_of_inertia * global_vector_angular_velocity / global_angular_velocity).norm();
        }

      std::vector<std::string> names;
      std::vector<std::string> units;

      names.push_back("Angular momentum");
      names.push_back("Moment of inertia");
      names.push_back("Angular velocity");

      if (this->convert_output_to_years() == true)
        {
          units.push_back("kg*m^2/year");
          units.push_back("kg*m^2");
          units.push_back("1/year");
          global_angular_momentum *= year_in_seconds;
          global_angular_velocity *= year_in_seconds;
        }
      else
        {
          units.push_back("kg*m^2/s");
          units.push_back("kg*m^2");
          units.push_back("1/s");
        }

      add_scientific_column(names[0] + " (" + units[0] +")", global_angular_momentum, statistics);

      if (dim == 2 || !output_full_tensor)
        add_scientific_column(names[1] + " (" + units[1] +")", global_moment_of_inertia, statistics);
      else
        {
          add_scientific_column(names[1] + "_xx (" + units[1] +")", global_tensor_moment_of_inertia[0][0], statistics);
          add_scientific_column(names[1] + "_yy (" + units[1] +")", global_tensor_moment_of_inertia[1][1], statistics);
          add_scientific_column(names[1] + "_zz (" + units[1] +")", global_tensor_moment_of_inertia[2][2], statistics);
          add_scientific_column(names[1] + "_xy (" + units[1] +")", global_tensor_moment_of_inertia[0][1], statistics);
          add_scientific_column(names[1] + "_xz (" + units[1] +")", global_tensor_moment_of_inertia[0][2], statistics);
          add_scientific_column(names[1] + "_yz (" + units[1] +")", global_tensor_moment_of_inertia[1][2], statistics);
        }

      add_scientific_column(names[2] + " (" + units[2] +")", global_angular_velocity, statistics);

      std::ostringstream output;
      output.precision(3);

      output << global_angular_momentum << " " << units[0] << ", "
             << global_moment_of_inertia << " " << units[1] << ", "
             << global_angular_velocity << " " << units[2];

      return std::pair<std::string, std::string> (names[0]+ ", " + names[1] + ", " + names[2] + ":",
                                                  output.str());
    }



    template <int dim>
    void
    RotationStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Rotation statistics");
        {
          prm.declare_entry("Use constant density of one","false",
                            Patterns::Bool(),
                            "Whether to use a constant density of one for the computation of the "
                            "angular momentum and moment of inertia. This is an approximation "
                            "that assumes that the 'volumetric' rotation is equal to the 'mass' "
                            "rotation. If this parameter is true this postprocessor computes "
                            "'net rotation' instead of 'angular momentum'.");
          prm.declare_entry("Output full moment of inertia tensor","false",
                            Patterns::Bool(),
                            "Whether to write the full moment of inertia tensor into the "
                            "statistics output instead of its norm for the current rotation "
                            "axis. This is a second-order symmetric tensor with "
                            "6 components in 3D. In 2D this option has no effect, because "
                            "the rotation axis is fixed and thus it is always a scalar.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    RotationStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Rotation statistics");
        {
          use_constant_density = prm.get_bool("Use constant density of one");
          output_full_tensor = prm.get_bool("Output full moment of inertia tensor");
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
    ASPECT_REGISTER_POSTPROCESSOR(RotationStatistics,
                                  "rotation statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "rotational velocity of the model (i.e. integrated "
                                  "net rotation and angular momentum.")
  }
}
