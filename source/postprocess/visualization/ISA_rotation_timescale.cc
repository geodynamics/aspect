// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/ISA_rotation_timescale.h>

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
      ISARotationTimescale<dim>::
      ISARotationTimescale ()
        :
        CellDataVectorCreator<dim>("s")
      {}



      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      ISARotationTimescale<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>> return_value("ISA_rotation_timescale",
                                                                              std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));

        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values(this->get_mapping(), this->get_fe(),
                                quadrature_formula,
                                update_values | update_gradients | update_quadrature_points);

        // Set up material models
        MaterialModel::MaterialModelInputs<dim> in(n_q_points,
                                                   this->n_compositional_fields());

        // Loop over cells and calculate tauISA in each one
        // Note that we start after timestep 0 because we need the strain rate,
        // which doesn't exist during the initial step
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned() && this->get_timestep_number() > 0)
            {

              // Fill the material model objects for the cell (for strain rate)
              fe_values.reinit(cell);
              in.reinit(fe_values, cell, this->introspection(),
                        this->get_solution());

              // Calculate eigenvalues of strain rate and take maximum (absolute value)
              // to get tauISA, the timescale for grain rotation toward the infinite strain axis
              const SymmetricTensor<2, dim> strain_rate = in.strain_rate[0];
              const std::array<double, dim> strain_rate_eigenvalues = eigenvalues(
                                                                        strain_rate);
              const double lambda1 = std::max(std::abs(strain_rate_eigenvalues[0]),
                                              std::abs(strain_rate_eigenvalues[dim-1]));
              const double tauISA = 1.0 / lambda1;

              (*return_value.second)(cell->active_cell_index()) = tauISA;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ISARotationTimescale,
                                                  "ISA rotation timescale",
                                                  "A visualization output object that generates output "
                                                  "showing the timescale for the rotation of grains "
                                                  "toward the infinite strain axis. Kaminski and Ribe "
                                                  "(see \\cite{Kaminski2002}) call this quantity "
                                                  "$\\tau_\\text{ISA}$ and define it as "
                                                  "$\\tau_\\text{ISA} \\approx \\frac{1}{\\dot{\\epsilon}}$ "
                                                  "where $\\dot{\\epsilon}$ is the largest eigenvalue "
                                                  "of the strain rate tensor. It can be used, "
                                                  "along with the grain lag angle $\\Theta$, "
                                                  "to calculate the grain orientation lag parameter."
                                                  "\n\n"
                                                  "Physical units: $\\text{s}$.")
    }
  }
}
