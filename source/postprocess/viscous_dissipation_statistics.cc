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
/*  $Id$  */


#include <aspect/postprocess/viscous_dissipation_statistics.h>
#include <aspect/simulator_access.h>
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
    ViscousDissipationStatistics<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim> quadrature_formula (this->get_fe()
                                            .base_element(this->introspection().base_elements.velocities)
                                            .degree+1);

      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_gradients |
                               update_JxW_values);

      double local_dissipation_integral = 0;

      // the values of the compositional fields are stored as blockvectors for each field
      // we have to extract them in this structure
      std::vector<std::vector<double> > prelim_composition_values (this->n_compositional_fields(),
                                                                   std::vector<double> (n_q_points));

      typename MaterialModel::Interface<dim>::MaterialModelInputs in(n_q_points,
                                                                     this->n_compositional_fields());
      typename MaterialModel::Interface<dim>::MaterialModelOutputs out(n_q_points,
                                                                       this->n_compositional_fields());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // retrieve the input for the material model
            fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                      in.pressure);
            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         in.temperature);
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values
              (this->get_solution(),prelim_composition_values[c]);

            for (unsigned int i=0; i<n_q_points; ++i)
              {
                for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                  in.composition[i][c] = prelim_composition_values[c][i];
              }

            fe_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(),
                in.strain_rate);

            // get the viscosity from the material model
            this->get_material_model().evaluate(in, out);

            // calculate the local viscous dissipation integral
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double div_v = trace(in.strain_rate[q]);
                local_dissipation_integral += ( - in.pressure[q] * div_v
                                                + 2.0 * out.viscosities[q] * in.strain_rate[q] * in.strain_rate[q]
                                                - (2.0 * out.viscosities[q] / 3.0) * div_v * div_v)
                                              * fe_values.JxW(q);
              }
          }

      // compute the viscous dissipation of the whole domain
      const double viscous_dissipation
        = Utilities::MPI::sum (local_dissipation_integral, this->get_mpi_communicator());

      if (this->convert_output_to_years() == true)
        {
          // fill statistics file
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.add_value ("Total viscous dissipation (J/yr)", viscous_dissipation * year_in_seconds);
          statistics.set_precision ("Total viscous dissipation (J/yr)", 8);
          statistics.set_scientific ("Total viscous dissipation (J/yr)", true);
        }
      else
        {
          // fill statistics file
          // make sure that the columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.add_value ("Total viscous dissipation (W)", viscous_dissipation);
          statistics.set_precision ("Total viscous dissipation (W)", 8);
          statistics.set_scientific ("Total viscous dissipation (W)", true);
        }

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << viscous_dissipation *year_in_seconds
               << " J/yr";
      else
        output << viscous_dissipation
               << " W";

      return std::pair<std::string, std::string> ("Total viscous dissipation:",
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
                                  "A postprocessor that computes the viscous dissipation"
                                  "for the whole domain as: "
                                  "$\\frac{1}{2} \\int_{V} \\sigma : \\dot{\\epsilon}dV$ "
                                  "= $\\int_{V} (-p\\nabla \\cdot u+2\\mu\\dot{\\epsilon}:\\dot{\\epsilon} "
                                  "- \\frac{2\\mu}{3}(\\nabla\\cdot u)^{2}) dV$. "
                                 )
  }
}
