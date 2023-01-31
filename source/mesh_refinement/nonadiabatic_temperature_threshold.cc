/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/nonadiabatic_temperature_threshold.h>
#include <aspect/adiabatic_conditions/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    NonadiabaticTemperatureThreshold<dim>::tag_additional_cells () const
    {
      // tag_additional_cells is executed before the equations are solved
      // for the very first time. If we do not have the finite element, we
      // do not have the temperature and just do nothing in this plugin.
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;

      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.temperature).get_unit_support_points());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values);

      std::vector<double> temperature_values (quadrature.size());
      const unsigned int n_dofs_per_cell = this->get_fe().base_element(this->introspection().base_elements.temperature).dofs_per_cell;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool refine = false;

            fe_values.reinit(cell);
            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         temperature_values);

            // if the nonadiabatic temperature exceeds the threshold, cell is marked for refinement
            for (unsigned int j=0; j<n_dofs_per_cell; ++j)
              {
                const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(fe_values.quadrature_point(j));

                double nonadiabatic_temperature = 0;
                if (temperature_anomaly_type == absolute_value)
                  nonadiabatic_temperature = std::abs(temperature_values[j] - adiabatic_temperature);
                else if (temperature_anomaly_type == positive_only)
                  nonadiabatic_temperature = temperature_values[j] - adiabatic_temperature;
                else if (temperature_anomaly_type == negative_only)
                  nonadiabatic_temperature = adiabatic_temperature - temperature_values[j];
                else
                  AssertThrow (false, ExcNotImplemented());

                if (nonadiabatic_temperature > threshold)
                  {
                    refine = true;
                    break;
                  }
              }

            if (refine)
              {
                cell->clear_coarsen_flag ();
                cell->set_refine_flag ();
              }
          }
    }

    template <int dim>
    void
    NonadiabaticTemperatureThreshold<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Nonadiabatic temperature threshold");
        {
          prm.declare_entry ("Threshold",
                             "100",
                             Patterns::Double (0.),
                             "A threshold that the nonadiabatic temperature "
                             "will be evaluated against. "
                             "Units: \\si{\\kelvin}");
          prm.declare_entry ("Temperature anomaly type",
                             "absolute value",
                             Patterns::Selection ("negative only|positive only|absolute value"),
                             "What type of temperature anomaly should be considered when "
                             "evaluating against the threshold: Only negative anomalies "
                             "(negative only), only positive anomalies (positive only) "
                             "or the absolute value of the nonadiabatic temperature.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    NonadiabaticTemperatureThreshold<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Nonadiabatic temperature threshold");
        {
          threshold = prm.get_double("Threshold");

          if (prm.get ("Temperature anomaly type") == "negative only")
            temperature_anomaly_type = negative_only;
          else if (prm.get ("Temperature anomaly type") == "positive only")
            temperature_anomaly_type = positive_only;
          else if (prm.get ("Temperature anomaly type") == "absolute value")
            temperature_anomaly_type = absolute_value;
          else
            AssertThrow (false, ExcNotImplemented());
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
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(NonadiabaticTemperatureThreshold,
                                              "nonadiabatic temperature threshold",
                                              "A mesh refinement criterion that computes refinement "
                                              "indicators from the temperature difference between the "
                                              "actual temperature and the adiabatic conditions (the "
                                              "nonadiabatic temperature). If the temperature anomaly "
                                              "exceeds the threshold given in the input file, the cell "
                                              "is marked for refinement.")
  }
}
