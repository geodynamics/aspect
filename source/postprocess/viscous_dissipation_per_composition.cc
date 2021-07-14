/*
  Copyright (C) 2011-2021 by the authors of the ASPECT code.

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


#include <aspect/postprocess/viscous_dissipation_per_composition.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ViscousDissipationPerComposition<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the compositional element alone.
      // be defensive about determining that a compositional field actually exists
      AssertThrow(this->introspection().base_elements.compositional_fields != numbers::invalid_unsigned_int,
                  ExcMessage("This postprocessor cannot be used without compositional fields."));
      const QGauss<dim> quadrature_formula(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree + 1);

      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_gradients |
                               update_JxW_values);

      std::vector<double> local_dissipation_integral (this->n_compositional_fields());

      std::vector<double> compositional_values(n_q_points);

      // the values of the compositional fields are stored as blockvectors for each field
      // we have to extract them in this structure
      std::vector<std::vector<double> > prelim_composition_values (this->n_compositional_fields(),
                                                                   std::vector<double> (n_q_points));

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points,
                                                                     this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_q_points,
                                                                       this->n_compositional_fields());

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            in.reinit(fe_values, cell, this->introspection(), this->get_solution());

            this->get_material_model().fill_additional_material_model_inputs(in, this->get_solution(), fe_values, this->introspection());

            this->get_material_model().evaluate(in, out);

            for (unsigned int c = 0; c<this->n_compositional_fields(); ++c)
              {
                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    if (in.composition[q][c] >= 0.5)
                      {
                        // Viscous dissipation D:
                        // D = 1/2 * volume_integral(deviatoric_stress*deviatoric_strain_rate)
                        //   = 1/2 * volume_integral(2*eta*deviatoric_strain_rate*deviatoric_strain_rate)
                        //   = volume_integral(eta*deviatoric_strain_rate*deviatoric_strain_rate)
                        const SymmetricTensor<2, dim> deviatoric_strain_rate =
                          (this->get_material_model().is_compressible()
                           ? in.strain_rate[q] - 1. / 3. * trace(in.strain_rate[q]) * unit_symmetric_tensor<dim>()
                           : in.strain_rate[q]);
                        local_dissipation_integral[c] += (out.viscosities[q] * deviatoric_strain_rate) *
                                                         deviatoric_strain_rate * fe_values.JxW(q);
                      }
                  }
              }
          }

      std::vector<double> viscous_dissipation (local_dissipation_integral.size());
      Utilities::MPI::sum (local_dissipation_integral, this->get_mpi_communicator(), viscous_dissipation);

      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          if (this->convert_output_to_years() == true)
            {
              statistics.add_value("Viscous dissipation (J/year/m) for composition " + this->introspection().name_for_compositional_index(c),
                                   viscous_dissipation[c] * year_in_seconds);

              // also make sure that the other columns filled by the this object
              // all show up with sufficient accuracy and in scientific notation
              {
                const std::string columns[] = { "Viscous dissipation (J/year/m) for composition " + this->introspection().name_for_compositional_index(c)};
                for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
                  {
                    statistics.set_precision (columns[i], 8);
                    statistics.set_scientific (columns[i], true);
                  }
              }
            }
          else
            {
              statistics.add_value ("Viscous dissipation (W/m) for composition " + this->introspection().name_for_compositional_index(c),
                                    viscous_dissipation[c]);

              // also make sure that the other columns filled by the this object
              // all show up with sufficient accuracy and in scientific notation
              {
                const std::string columns[] = { "Viscous dissipation (W/m) for composition " + this->introspection().name_for_compositional_index(c)};
                for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
                  {
                    statistics.set_precision (columns[i], 8);
                    statistics.set_scientific (columns[i], true);
                  }
              }
            }
        }

      std::ostringstream output;
      output.precision(4);
      for (unsigned int c = 0; c < this->n_compositional_fields(); ++c)
        {
          if (this->convert_output_to_years() == true)
            output << viscous_dissipation[c] * year_in_seconds
                   << " J/year/m";
          else
            output << viscous_dissipation[c]
                   << " W/m";

          if (c + 1 != this->n_compositional_fields())
            output << " // ";
        }

      return std::pair<std::string, std::string> ("Compositions viscous dissipation:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ViscousDissipationPerComposition,
                                  "viscous dissipation per compositional field",
                                  "A postprocessor that outputs the viscous rate of dissipation of "
                                  "energy for each compositional field (where the field has a value "
                                  "of 0.5 or more). "
                                  "The viscous dissipation is computed as: "
                                  "$\\int_{V}\\left(\\sigma' \\dot{\\epsilon}' \\right)$, "
                                  "where $\\sigma'$  is the deviatoric stress and $\\dot{\\epsilon}'$ "
                                  "the deviatoric strain rate.")
  }
}
