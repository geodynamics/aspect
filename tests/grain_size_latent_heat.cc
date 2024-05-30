/*
  Copyright (C) 2014 - 2023 by the authors of the ASPECT code.

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
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except that the density decays linearly with the
     * temperature and the viscosity, which depends on the temperature,
     * pressure, strain rate and grain size.
     *
     * The grain size evolves in time, dependent on strain rate, temperature,
     * creep regime, and phase transitions.
     *
     * The model is considered compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class GrainSizeLatentHeat : public MaterialModel::GrainSize<dim>
    {
      public:
        virtual bool is_compressible () const override
        {
          return false;
        }

        virtual void evaluate(const typename MaterialModel::Interface<dim>::MaterialModelInputs &in,
                              typename MaterialModel::Interface<dim>::MaterialModelOutputs &out) const override
        {
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
                  const double own_enthalpy = this->material_lookup[0]->enthalpy(temperatures[q],pressures[q]);
                  for (unsigned int p=0; p<n_q_points; ++p)
                    {
                      double enthalpy_p,enthalpy_T;
                      if (std::fabs(temperatures[q] - temperatures[p]) > 1e-12 * temperatures[q])
                        {
                          enthalpy_p = this->material_lookup[0]->enthalpy(temperatures[p],pressures[q]);
                          const double point_contribution = (own_enthalpy-enthalpy_p)/(temperatures[q]-temperatures[p]);
                          dHdT += point_contribution;
                          T_points++;
                        }
                      if (std::fabs(pressures[q] - pressures[p]) > 1)
                        {
                          enthalpy_T = this->material_lookup[0]->enthalpy(temperatures[q],pressures[p]);
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
              // convert the grain size from log to normal
              std::vector<double> composition (in.composition[i]);
              if (this->advect_log_grainsize)
                this->convert_log_grain_size(composition);
              else
                for (unsigned int c=0; c<composition.size(); ++c)
                  composition[c] = std::max(this->min_grain_size,composition[c]);

              // Set up an integer that tells us which phase transition has been crossed inside of the cell.
              int crossed_transition(-1);

              // Figure out if the material in the current cell underwent a phase change.
              // To do so, check if a grain has moved further than the distance from the phase transition and
              // if the velocity is in the direction of the phase change. After the check 'crossed_transition' will
              // be -1 if we crossed no transition, or the index of the phase transition, if we crossed it.
              for (unsigned int phase=0; phase<this->n_phase_transitions; ++phase)
                {
                  const Tensor<1,dim> vertical_direction = this->get_gravity_model().gravity_vector(in.position[i])
                                                           /this->get_gravity_model().gravity_vector(in.position[i]).norm();
                  const double timestep = this->simulator_is_past_initialization()
                                          ?
                                          this->get_timestep()
                                          :
                                          0.0;

                  // Both distances are positive when they are downward from the transition (since gravity points down)
                  const double distance_from_transition = this->get_geometry_model().depth(in.position[i]) - this->phase_function.get_transition_depth(phase);
                  const double distance_moved = in.velocity[i] * vertical_direction * timestep;

                  // If we are close to the phase boundary (closer than the distance a grain has moved
                  // within one time step) and the velocity points away from the phase transition,
                  // then the material has crossed the transition.
                  // To make sure we actually reset the grain size of all the material passing through
                  // the transition, we take 110% of the distance a grain has moved for the check.
                  if (std::abs(distance_moved) * 1.1 > std::abs(distance_from_transition)
                      &&
                      distance_moved * distance_from_transition >= 0)
                    crossed_transition = phase;
                }

              if (in.requests_property(MaterialProperties::viscosity))
                {
                  double effective_viscosity;
                  double disl_viscosity = std::numeric_limits<double>::max();
                  Assert(std::isfinite(in.strain_rate[i].norm()),
                         ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                    "not filled by the caller."));
                  const SymmetricTensor<2,dim> shear_strain_rate = in.strain_rate[i] - 1./dim * trace(in.strain_rate[i]) * unit_symmetric_tensor<dim>();
                  const double second_strain_rate_invariant = std::sqrt(std::max(-second_invariant(shear_strain_rate), 0.));

                  const double adiabatic_temperature = this->get_adiabatic_conditions().is_initialized()
                                                       ?
                                                       this->get_adiabatic_conditions().temperature(in.position[i])
                                                       :
                                                       in.temperature[i];
                  const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                                    ?
                                                    this->get_adiabatic_conditions().pressure(in.position[i])
                                                    :
                                                    in.pressure[i];

                  const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");

                  const double diff_viscosity = this->diffusion_viscosity(in.temperature[i],
                                                                          adiabatic_temperature,
                                                                          adiabatic_pressure,
                                                                          composition[grain_size_index],
                                                                          second_strain_rate_invariant,
                                                                          in.position[i]);

                  if (std::abs(second_strain_rate_invariant) > 1e-30)
                    {
                      disl_viscosity = this->dislocation_viscosity(in.temperature[i], adiabatic_temperature, adiabatic_pressure, in.strain_rate[i], in.position[i],diff_viscosity);
                      effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
                    }
                  else
                    effective_viscosity = diff_viscosity;

                  out.viscosities[i] = std::min(std::max(this->min_eta,effective_viscosity),this->max_eta);
                }

              out.densities[i] = this->density(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);

              if (this->get_adiabatic_conditions().is_initialized())
                {
                  if ((in.current_cell.state() == IteratorState::valid)
                      && (std::fabs(dHdp) > std::numeric_limits<double>::epsilon())
                      && (std::fabs(dHdT) > std::numeric_limits<double>::epsilon()))
                    {
                      out.thermal_expansion_coefficients[i] = (1 - 3515.6 * dHdp) / in.temperature[i];
                      out.specific_heat[i] = dHdT;
                    }
                  else
                    {
                      out.thermal_expansion_coefficients[i] = this->thermal_expansion_coefficient(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
                      out.specific_heat[i] = this->specific_heat(in.temperature[i], in.pressure[i], composition, in.position[i]);
                    }
                }
              else
                {
                  out.thermal_expansion_coefficients[i] = (1 - 3515.6 * this->material_lookup[0]->dHdp(in.temperature[i],in.pressure[i])) / in.temperature[i];
                  out.specific_heat[i] = this->material_lookup[0]->dHdT(in.temperature[i],in.pressure[i]);
                }

              out.thermal_conductivities[i] = this->k_value;
              out.compressibilities[i] = this->compressibility(in.temperature[i], in.pressure[i], composition, in.position[i]);

              // TODO: make this more general for not just olivine grains
              if (in.requests_property(MaterialProperties::reaction_terms))
                for (unsigned int c=0; c<composition.size(); ++c)
                  {
                    if (this->introspection().name_for_compositional_index(c) == "olivine_grain_size")
                      {
                        out.reaction_terms[i][c] = this->grain_size_change(in.temperature[i], in.pressure[i], composition,
                                                                           in.strain_rate[i], in.velocity[i], in.position[i], c, crossed_transition);
                        if (this->advect_log_grainsize)
                          out.reaction_terms[i][c] = - out.reaction_terms[i][c] / composition[c];
                      }
                    else
                      out.reaction_terms[i][c] = 0.0;
                  }
            }
        }
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
