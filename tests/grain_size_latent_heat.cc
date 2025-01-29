/*
  Copyright (C) 2014 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/grain_size.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

#include <deal.II/base/parameter_handler.h>

#include <iostream>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that behaves in the same way as the grain size model, but is modified to
     * resemble the latent heat benchmark. Due to the nature of the benchmark the model needs to be
     * incompressible despite using a material table. It assumes a constant density for the calculation of
     * the latent heat.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class GrainSizeLatentHeat : public MaterialModel::Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        bool is_compressible () const override
        {
          return false;
        }

        void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                      typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const override
        {
          base_model->evaluate(in, out);

          double dHdT = 0.0;
          double dHdp = 0.0;

          std::vector<double> compositional_fields(this->n_compositional_fields(), 0.);
          compositional_fields[0] = 1.0;

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

              unsigned int T_points(0),p_points(0);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double own_enthalpy = base_model->enthalpy(temperatures[q],pressures[q],compositional_fields,Point<dim>());
                  for (unsigned int p=0; p<n_q_points; ++p)
                    {
                      double enthalpy_p,enthalpy_T;
                      if (std::fabs(temperatures[q] - temperatures[p]) > 1e-12 * temperatures[q])
                        {
                          enthalpy_p = base_model->enthalpy(temperatures[p],pressures[q],compositional_fields,Point<dim>());
                          const double point_contribution = (own_enthalpy-enthalpy_p)/(temperatures[q]-temperatures[p]);
                          dHdT += point_contribution;
                          T_points++;
                        }
                      if (std::fabs(pressures[q] - pressures[p]) > 1)
                        {
                          enthalpy_T = base_model->enthalpy(temperatures[q],pressures[p],compositional_fields,Point<dim>());
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
              const double approximate_density = 3515.6;

              if (this->get_adiabatic_conditions().is_initialized())
                {
                  if ((in.current_cell.state() == IteratorState::valid)
                      && (std::fabs(dHdp) > std::numeric_limits<double>::epsilon())
                      && (std::fabs(dHdT) > std::numeric_limits<double>::epsilon()))
                    {
                      out.thermal_expansion_coefficients[i] = (1 - approximate_density * dHdp) / in.temperature[i];
                    }
                }
              else
                {
                  // Estimate the thermal expansivity by approximating dHdp at constant
                  // temperature from the enthalpy lookup.
                  // The data table has a pressure increment of 8e8 Pa.
                  const double delta_pressure = 8e8;

                  // compositional fields and position are not used for this test in the base model
                  const double h = base_model->enthalpy(in.temperature[i],in.pressure[i],compositional_fields,Point<dim>());
                  const double dh = base_model->enthalpy(in.temperature[i],in.pressure[i] + delta_pressure,compositional_fields,Point<dim>());
                  dHdp = (dh - h) / delta_pressure;

                  out.thermal_expansion_coefficients[i] = (1 - approximate_density * dHdp) / in.temperature[i];
                }
            }
        }

        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          MaterialModel::GrainSize<dim>::declare_parameters(prm);
        }

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override
        {
          base_model = std::make_unique<MaterialModel::GrainSize<dim>>();
          base_model->initialize_simulator(this->get_simulator());
          base_model->parse_parameters(prm);
          base_model->initialize();
        }

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override
        {
          base_model->create_additional_named_outputs(out);
        }

      private:
        std::unique_ptr<MaterialModel::GrainSize<dim>> base_model;
    };
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(GrainSizeLatentHeat,
                                   "grain size latent heat",
                                   "A material model that behaves in the same way as "
                                   "the grain size model, but is modified to "
                                   "resemble the latent heat benchmark. Due to the "
                                   "nature of the benchmark the model needs to be "
                                   "incompressible despite a material table and "
                                   "use a constant density for the calculation of "
                                   "the latent heat.")
  }
}
