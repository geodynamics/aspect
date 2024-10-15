/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/grain_lag_angle.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/symmetric_tensor.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      GrainLagAngle<dim>::
      GrainLagAngle ()
        :
        CellDataVectorCreator<dim>("radian")
      {}



      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      GrainLagAngle<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>> return_value("grain_lag_angle",
                                                                              std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));

        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size(); // this is 1 for QMidpoint

        FEValues<dim> fe_values(this->get_mapping(), this->get_fe(),
                                quadrature_formula,
                                update_values | update_gradients | update_quadrature_points);

        // Set up material models
        MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                   this->n_compositional_fields());

        // Loop over cells and calculate theta in each one
        // Note that we start after timestep 0 because we need the strain rate,
        // which doesn't exist during the initial step
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned() && this->get_timestep_number() > 0)
            {
              // Fill the material model objects for the cell (for strain rate)
              fe_values.reinit(cell);
              in.reinit(fe_values, cell, this->introspection(),
                        this->get_solution());
              // Also get velocity gradients
              std::vector<Tensor<2, dim>> velocity_gradient(n_q_points,
                                                             Tensor<2, dim>());
              fe_values[this->introspection().extractors.velocities].get_function_gradients(
                this->get_solution(), velocity_gradient);

              // Calculate eigenvalues of strain rate and take maximum (absolute value)
              // to get tauISA, the timescale for grain rotation toward the infinite strain axis
              const SymmetricTensor<2, dim> strain_rate = in.strain_rate[0];
              const std::array<double, dim> strain_rate_eigenvalues = eigenvalues(
                                                                        strain_rate);

              const double lambda1 = std::max(std::abs(strain_rate_eigenvalues[0]),
                                              std::abs(strain_rate_eigenvalues[dim-1]));
              const double tauISA = 1.0 / lambda1;

              // Next, we want the ~infinite strain axis, so we start by
              // calculating the deformation gradient tensor in the limit of "infinite"
              // time using the velocity gradient tensor. The deformation gradient tensor (F) obeys
              // dF/dt = LF
              // where L is the velocity gradient tensor, and F(t=0) = I. The solution to this is
              // F = exp(Lt)
              // which can be written as an infinite sum:
              // F = I + Lt + L^2*t^2/2! + L^3*t^3/3! + ...
              // and then the "left stretch" tensor is U = F * F^T
              // We truncate the series sum at the L^4 term. Note that L is sometimes
              // taken to be the strain rate tensor (ie the velocity gradient without
              // the rotational component).
              // The eigenvectors of U are the directions of the major axes of the finite
              // strain ellipsoid, so for the *infinite* strain axis, we want the eigenvector
              // corresponding to the largest eigenvalue of U in the limit as t->infinity.
              // As per Kaminski and Ribe (2002), tmax = 75*tauISA is infinite enough for
              // most applications.
              const unsigned int maxorder = 4; // truncate at L^4
              Tensor<2, dim> F = unit_symmetric_tensor<dim, double>(); // initialize as identity (first term of sum)
              const double tmax = 75 * tauISA; // "infinite" time

              for (unsigned int i = 1; i <= maxorder; ++i) // loop and calculate terms of sum
                {
                  Tensor<2, dim> multiplier = unit_symmetric_tensor<dim, double>();  // to collect powers of L
                  double factor = 1.0;  // to collect powers of tmax and factorials
                  for (unsigned int j = 0; j < i; ++j)
                    {
                      multiplier = multiplier
                                   * Tensor<2, dim>(velocity_gradient[0]);
                      factor *= tmax / (j + 1);
                    }
                  multiplier = factor * multiplier;
                  F += multiplier;  // sum terms
                }

              // The "left stretch" tensor is computed from the deformation gradient tensor
              // Note that Kaminski and Ribe incorrectly wrote U=Ft*F (typo, I guess).
              const Tensor<2, dim> Ft = transpose(F);
              const SymmetricTensor<2, dim> U = SymmetricTensor<2, dim>(F * Ft);
              // The ~infinite strain axis is the first eigenvector of the "left stretch" tensor.
              const Tensor<1, dim> ehat = eigenvectors(U)[0].second;
              // Calculate theta by dotting the ISA (ehat) and the velocity (scaled to a unit vector)
              // theta is limited to be less than pi/2 (stay in the first quadrant)
              const double umag = in.velocity[0].norm();
              double theta_val = 0;
              if (umag > 100.0 * std::numeric_limits<double>::min())
                {
                  double dotprod = 0;
                  for (unsigned int iq = 0; iq < dim; ++iq)
                    {
                      dotprod = dotprod + in.velocity[0][iq]*ehat[iq];
                    }
                  theta_val = std::acos(std::abs(dotprod / umag));
                  if (theta_val > numbers::PI/2)
                    theta_val = std::fmod(theta_val, numbers::PI/2);
                }
              (*return_value.second)(cell->active_cell_index()) = theta_val;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(GrainLagAngle,
                                                  "grain lag angle",
                                                  "A visualization output object that generates output "
                                                  "showing the angle between the ~infinite strain axis "
                                                  "and the flow velocity. Kaminski and Ribe "
                                                  "(see \\cite{Kaminski2002}) call this quantity $\\Theta$ and "
                                                  "define it as "
                                                  "$\\Theta = \\cos^{-1}(\\hat{u}\\cdot\\hat{e})$ "
                                                  " where $\\hat{u}=\\vec{u}/|{u}|$, $\\vec{u}$ "
                                                  "is the local flow velocity, and "
                                                  "$\\hat{e}$ is the local infinite strain axis, "
                                                  "which we calculate as the first eigenvector of "
                                                  "the 'left stretch' tensor. "
                                                  "$\\Theta$ can be used to calculate the grain "
                                                  "orientation lag parameter."
                                                  "\n\n"
                                                  "Physical units: \\si{\\radian}.")
    }
  }
}
