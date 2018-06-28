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


#include <aspect/postprocess/visualization/theta_tauisa.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/symmetric_tensor.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template<int dim>
      std::pair<std::string, Vector<float> *> Theta<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *> return_value("theta",
                                                             new Vector<float>(this->get_triangulation().n_active_cells()));

        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size(); // this is 1 for QMidpoint

        FEValues<dim> fe_values(this->get_mapping(), this->get_fe(),
                                quadrature_formula,
                                update_values | update_gradients | update_quadrature_points);

        // Set up material models
        MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                     this->n_compositional_fields());

        // Set up cell iterator for looping
        typename DoFHandler<dim>::active_cell_iterator cell =
          this->get_dof_handler().begin_active(), endc =
            this->get_dof_handler().end();

        // Loop over cells and calculate theta in each one
        // Note that we start after timestep 0 because we need the strain rate,
        // which doesn't exist during the initial step
        unsigned int cell_index = 0;
        for (; cell != endc; ++cell, ++cell_index)
          {
            if (cell->is_locally_owned() && this->get_timestep_number() > 0)
              {

                // Fill the material model objects for the cell (for strain rate)
                fe_values.reinit(cell);
                in.reinit(fe_values, cell, this->introspection(),
                          this->get_solution(), true);
                // Also get velocities and velocity gradients
                std::vector<Tensor<1, dim> > velocities(n_q_points,
                                                        Tensor<1, dim>());
                std::vector<Tensor<2, dim> > velocity_gradient(n_q_points,
                                                               Tensor<2, dim>());
                fe_values[this->introspection().extractors.velocities].get_function_values(
                  this->get_solution(), velocities);
                fe_values[this->introspection().extractors.velocities].get_function_gradients(
                  this->get_solution(), velocity_gradient);

                // Calculate eigenvalues of strain rate and take maximum (absolute value)
                // to get tauISA, the timescale for grain rotation toward the infinite strain axis
                const SymmetricTensor<2, dim> strain_rate = in.strain_rate[0];
                std::array<double, dim> strain_rate_eigenvalues = eigenvalues(
                                                                    strain_rate);
                double lamda1 = std::max(std::abs(strain_rate_eigenvalues[0]),
                                         std::abs(strain_rate_eigenvalues[dim]));
                double tauISA = 1.0 / lamda1;

                // Get the deformation gradient tensor in the limit of "infinite" time
                // from the velocity gradient tensor. As per Kaminski and Ribe (2002),
                // tmax = 75*tauISA is infinite enough for most applications.
                unsigned int maxorder = 4;
                Tensor<2, dim> F = Tensor<2, dim>(
                                     unit_symmetric_tensor<dim, double>());
                double tmax = 75 * tauISA;

                for (unsigned int i = 0; i < maxorder; ++i)
                  {
                    Tensor<2, dim> multiplier = Tensor<2, dim>(
                                                  unit_symmetric_tensor<dim, double>());
                    double factor = 1.0;
                    for (unsigned int j = 0; j < i; ++j)
                      {
                        multiplier = multiplier
                                     * Tensor<2, dim>(velocity_gradient[0]);
                        factor *= tmax / (j + 1);
                      }
                    multiplier = factor * multiplier;
                    F += multiplier;
                  }

                // The "left stretch" tensor is computed from the deformation gradient tensor
                // Note that Kaminski and Ribe incorrectly wrote U=Ft*F (typo, I guess).
                Tensor<2, dim> Ft = transpose(F);
                SymmetricTensor<2, dim> U = SymmetricTensor<2, dim>(F * Ft);
                // The ~infinite strain axis is the first eigenvector of the "left stretch" tensor.
                Tensor<1, dim> ehat = eigenvectors(U)[0].second;

                // Calculate theta by dotting the ISA (ehat) and the velocity (scaled to a unit vector)
                // theta is limited to be less than pi/2 (stay in the first quadrant)
                double umag = std::sqrt(velocities[0] * velocities[0]);
                double umagmin = 1.0e-40;  // to avoid problems with /0 in acos
                double theta_val = 0;
                if (umag > umagmin)
                  {
                    theta_val = std::acos(std::abs((velocities[0] * ehat) / umag));
                    if (theta_val > 1.57079632679)
                      theta_val = fmod(theta_val, 1.57079632679);
                  }
                (*return_value.second)(cell_index) = theta_val;
              }
          }
        return return_value;
      }

      template<int dim>
      std::pair<std::string, Vector<float> *> tauISA<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *> return_value("tauISA",
                                                             new Vector<float>(this->get_triangulation().n_active_cells()));

        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values(this->get_mapping(), this->get_fe(),
                                quadrature_formula,
                                update_values | update_gradients | update_quadrature_points);

        // Set up material models
        MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                   this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_q_points,
                                                     this->n_compositional_fields());

        // Set up cell iterator for looping
        typename DoFHandler<dim>::active_cell_iterator cell =
          this->get_dof_handler().begin_active(), endc =
            this->get_dof_handler().end();

        // Loop over cells and calculate tauISA in each one
        // Note that we start after timestep 0 because we need the strain rate,
        // which doesn't exist during the initial step
        unsigned int cell_index = 0;
        for (; cell != endc; ++cell, ++cell_index)
          {
            if (cell->is_locally_owned() && this->get_timestep_number() > 0)
              {

                // Fill the material model objects for the cell (for strain rate)
                fe_values.reinit(cell);
                in.reinit(fe_values, cell, this->introspection(),
                          this->get_solution(), true);

                // Calculate eigenvalues of strain rate and take maximum (absolute value)
                // to get tauISA, the timescale for grain rotation toward the infinite strain axis
                const SymmetricTensor<2, dim> strain_rate = in.strain_rate[0];
                std::array<double, dim> strain_rate_eigenvalues = eigenvalues(
                                                                    strain_rate);
                double lamda1 = std::max(std::abs(strain_rate_eigenvalues[0]),
                                         std::abs(strain_rate_eigenvalues[dim]));
                double tauISA = 1.0 / lamda1;

                (*return_value.second)(cell_index) = tauISA;
              }
          }

        return return_value;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Theta, "Theta",
                                                  "A visualization output object that generates output"
                                                  "showing the angle between the infinite strain axis "
                                                  "(approximated by a truncated series solution) "
                                                  "and the flow velocity.")

      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(tauISA, "tauISA",
                                                  "A visualization output object that generates output "
                                                  "showing the timescale for the rotation of grains "
                                                  "toward the infinite strain axis.")
    }
  }
}
