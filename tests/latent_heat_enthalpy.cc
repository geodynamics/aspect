/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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


#include <aspect/material_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <iostream>
#include <deal.II/base/quadrature_lib.h>
#include <memory>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class LatentHeatEnthalpy : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void
        initialize() override
        {
          const std::string datadirectory = Utilities::expand_ASPECT_SOURCE_DIR("$ASPECT_SOURCE_DIR/data/material-model/latent-heat-enthalpy-test/");
          const std::string material_file_names  = "testdata.txt";

          material_lookup = std::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>(datadirectory+material_file_names,
                            true,
                            this->get_mpi_communicator());
        }

        bool
        is_compressible () const override
        {
          return false;
        }

        void
        evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const override
        {
          const double reference_rho              = 3400;
          const double density_jump               = 115.6;
          const double eta                        = 8.44e21;
          const double k_value                    = 2.38;

          double dHdT = 0.0;
          double dHdp = 0.0;

          if (in.current_cell.state() == IteratorState::valid)
            {
              const QTrapezoid<dim> quadrature_formula;
              const unsigned int n_q_points = quadrature_formula.size();

              FEValues<dim> fe_values (this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula,
                                       update_values);

              std::vector<double> temperatures(n_q_points), pressures(n_q_points);
              std::vector<std::vector<double>> compositions (quadrature_formula.size(),std::vector<double> (this->n_compositional_fields()));
              std::vector<std::vector<double>> composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

              fe_values.reinit (in.current_cell);

              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature]
              .get_function_values (this->get_solution(), temperatures);
              fe_values[this->introspection().extractors.pressure]
              .get_function_values (this->get_solution(), pressures);

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]]
                .get_function_values(this->get_solution(),
                                     composition_values[c]);
              for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
                {
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    compositions[q][c] = composition_values[c][q];
                }

              unsigned int T_points(0),p_points(0);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double own_enthalpy = material_lookup->enthalpy(temperatures[q],pressures[q]);
                  for (unsigned int p=0; p<n_q_points; ++p)
                    {
                      double enthalpy_p,enthalpy_T;
                      if (std::fabs(temperatures[q] - temperatures[p]) > 100 * std::fabs(temperatures[q]) * std::numeric_limits<double>::epsilon())
                        {
                          enthalpy_p = material_lookup->enthalpy(temperatures[p],pressures[q]);
                          const double point_contribution = (own_enthalpy-enthalpy_p)/(temperatures[q]-temperatures[p]);
                          dHdT += point_contribution;
                          T_points++;
                        }
                      if (std::fabs(pressures[q] - pressures[p]) > 100 * std::fabs(pressures[q]) * std::numeric_limits<double>::epsilon())
                        {
                          enthalpy_T = material_lookup->enthalpy(temperatures[q],pressures[p]);
                          dHdp += (own_enthalpy-enthalpy_T)/(pressures[q]-pressures[p]);
                          p_points++;
                        }
                    }
                }

              if ((T_points > 0)
                  && (p_points > 0))
                {
                  dHdT /= T_points;
                  dHdp /= p_points;
                }
            }

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              if (in.requests_property(MaterialProperties::viscosity))
                out.viscosities[i] = eta;

              out.densities[i] = material_lookup->density(in.temperature[i],in.pressure[i]);

              if (this->get_adiabatic_conditions().is_initialized())
                {
                  if ((in.current_cell.state() == IteratorState::valid)
                      && (std::fabs(dHdp) > 100 * std::numeric_limits<double>::epsilon())
                      && (std::fabs(dHdT) > 100 * std::numeric_limits<double>::epsilon()))
                    {
                      out.thermal_expansion_coefficients[i] = (1 - (reference_rho + density_jump) * dHdp) / in.temperature[i];
                      out.specific_heat[i] = dHdT;
                    }
                  else
                    {
                      out.thermal_expansion_coefficients[i] = material_lookup->thermal_expansivity(in.temperature[i],in.pressure[i]);
                      out.specific_heat[i] = material_lookup->specific_heat(in.temperature[i],in.pressure[i]);
                    }
                }
              else
                {
                  out.thermal_expansion_coefficients[i] = (1 - (reference_rho + density_jump) * material_lookup->dHdp(in.temperature[i],in.pressure[i])) / in.temperature[i];
                  out.specific_heat[i] = material_lookup->dHdT(in.temperature[i],in.pressure[i]);
                }

              out.thermal_conductivities[i] = k_value;
            }
        }

      private:
        std::unique_ptr<MaterialModel::MaterialUtilities::Lookup::MaterialLookup> material_lookup;
    };
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(LatentHeatEnthalpy,
                                   "latent heat enthalpy",
                                   "A test model that implements a latent "
                                   "heat formulation based on a lookup table "
                                   "that is evaluated at the edges of the cell "
                                   "rather than at the quadrature points.")
  }
}
