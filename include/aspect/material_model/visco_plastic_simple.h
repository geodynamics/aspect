/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_visco_plastic_simple_h
#define _aspect_material_model_visco_plastic_simple_h

#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
#include <aspect/material_model/rheology/diffusion_creep.h>
#include <aspect/material_model/rheology/dislocation_creep.h>
#include <aspect/material_model/rheology/elasticity.h>
#include <aspect/material_model/rheology/drucker_prager.h>
#include <aspect/material_model/rheology/friction_models.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class ViscoPlasticSimple : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        bool is_compressible() const override;

        static
        void
        declare_parameters(ParameterHandler &prm);

        void
        parse_parameters(ParameterHandler &prm) override;

        void
        create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:
        std::pair<double, double>
        calculate_viscosity_and_plastic_strain_rate(const unsigned int chemical_field,
                                                    const double strain_rate,
                                                    const double pressure,
                                                    const double temperature,
                                                    const double plastic_strain) const;

        double 
        cellwise_average(const MaterialAveraging::AveragingOperation operation,
                         const std::vector<double> &values) const;

        double ref_strain_rate;

        double min_strain_rate;

        double min_viscosity;

        double max_viscosity;

        bool define_conductivities;

        std::vector<double> thermal_conductivities;

        std::vector<double> thermal_diffusivities;

        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;

        MaterialUtilities::CompositionalAveragingOperation pointwise_viscosity_averaging;

        MaterialAveraging::AveragingOperation cellwise_viscosity_averaging;

        bool allow_negative_pressure_in_plasticity;

        bool enable_plastic_strain_weakening;

        std::vector<double> start_plastic_strain_weakening_intervals;

        std::vector<double> end_plastic_strain_weakening_intervals;

        std::vector<double> cohesion_strain_weakening_factors;

        std::vector<double> friction_strain_weakening_factors;

        enum ViscousFlowLaw
        {
          diffusion,
          dislocation,
          composite
        } viscous_flow_law;

        double newton_iteration_tolerance;

        unsigned int max_newton_iteration;

        Rheology::DiffusionCreep<dim> diffusion_creep;

        Rheology::DislocationCreep<dim> dislocation_creep;

        Rheology::Elasticity<dim> elastic_rheology;

        Rheology::DruckerPrager<dim> drucker_prager_plasticity;

        Rheology::FrictionModels<dim> friction_models;
    };
  }
}

#endif
