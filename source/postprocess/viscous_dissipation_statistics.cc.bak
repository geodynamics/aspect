/*
  Copyright (C) 2011 - 2023-2021 by the authors of the ASPECT code.

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


#include <aspect/postprocess/viscous_dissipation_statistics.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ViscousDissipationStatistics<dim>::execute (TableHandler &statistics)
    {
      // Create a quadrature formula based on the velocity element.
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.velocities;

      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_gradients |
                               update_JxW_values);

      // Vector to store the local dissipation for all the fields and
      // for the whole domain.
      const unsigned int n_compositional_fields = this->n_compositional_fields();

      std::vector<double> local_dissipation_integral (n_compositional_fields+1);
      std::vector<double> compositional_values(n_q_points);

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points,
                                                                     n_compositional_fields);
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_q_points,
                                                                       n_compositional_fields);
      in.requested_properties = MaterialModel::MaterialProperties::viscosity;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            in.reinit(fe_values,
                      cell,
                      this->introspection(),
                      this->get_solution());

            this->get_material_model().fill_additional_material_model_inputs(in,
                                                                             this->get_solution(),
                                                                             fe_values,
                                                                             this->introspection());

            this->get_material_model().evaluate(in, out);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                // Viscous dissipation D in 3d:
                // D = 1/2 * volume_integral(deviatoric_stress*deviatoric_strain_rate)
                //   = 1/2 * volume_integral(2*eta*deviatoric_strain_rate*deviatoric_strain_rate)
                //   = volume_integral(eta*deviatoric_strain_rate*deviatoric_strain_rate)
                const SymmetricTensor<2, dim> deviatoric_strain_rate =
                  (this->get_material_model().is_compressible()
                   ? in.strain_rate[q] - 1. / 3. * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>()
                   : in.strain_rate[q]);
                const double local_dissipation = out.viscosities[q] * (deviatoric_strain_rate *
                                                                       deviatoric_strain_rate) * fe_values.JxW(q);

                // Dissipation over the whole domain
                local_dissipation_integral[n_compositional_fields] += local_dissipation;

                // Dissipation over each compositional field
                for (unsigned int c = 0; c<n_compositional_fields; ++c)
                  if (in.composition[q][c] >= 0.5)
                    local_dissipation_integral[c] += local_dissipation;
              }
          }

      std::vector<double> viscous_dissipation (local_dissipation_integral.size());
      Utilities::MPI::sum (local_dissipation_integral, this->get_mpi_communicator(), viscous_dissipation);

      const std::string unit = (dim == 3) ? "(W)" : "(W/m)";

      // Add the dissipation per compositional fields to the statistics output
      for (unsigned int c = 0; c < n_compositional_fields; ++c)
        {
          statistics.add_value ("Viscous dissipation " + unit + " for composition " + this->introspection().name_for_compositional_index(c),
                                viscous_dissipation[c]);

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          const std::string column = "Viscous dissipation " + unit + " for composition " + this->introspection().name_for_compositional_index(c);

          statistics.set_precision (column, 8);
          statistics.set_scientific (column, true);
        }

      // Add the dissipation over the whole domain to the statistics output
      statistics.add_value ("Viscous dissipation " + unit,
                            viscous_dissipation[n_compositional_fields]);

      const std::string column = "Viscous dissipation " + unit;
      statistics.set_precision (column, 8);
      statistics.set_scientific (column, true);

      std::ostringstream output;
      output.precision(4);

      for (unsigned int n = 0; n < n_compositional_fields + 1; ++n)
        {
          output << viscous_dissipation[n];

          if (n + 1 != n_compositional_fields+1)
            output << " // ";
        }

      return std::pair<std::string, std::string> ("Viscous dissipation per field and for whole domain " + unit + ":",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ViscousDissipationStatistics,
                                  "viscous dissipation statistics",
                                  "A postprocessor that outputs the viscous rate of dissipation of "
                                  "energy for each compositional field (where the field has a value "
                                  "of 0.5 or more) as well as over the whole domain. "
                                  "When all the fields represent lithologies and there is no background "
                                  "field, the sum of the individual field's dissipation should equal "
                                  "that over the whole domain. "
                                  "The viscous dissipation is computed as: "
                                  "$\\int_{V}\\left(\\sigma^\\prime \\dot{\\epsilon}^\\prime \\right)$, "
                                  "where $\\sigma^\\prime$  is the deviatoric stress and $\\dot{\\epsilon}^\\prime$ "
                                  "the deviatoric strain rate."
                                  "Note then when shear heating is included in the temperature equation, "
                                  "it is better to use the 'heating statistics' postprocessor.")
  }
}
