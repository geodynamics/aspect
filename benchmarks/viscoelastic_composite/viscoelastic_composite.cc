/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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


#include "viscoelastic_composite.h"
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>
#include <aspect/material_model/equation_of_state/interface.h>

#include <aspect/heating_model/shear_heating.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    ViscoelasticComposite<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // The following is the copy from the Simpler material model
      // because I can't think of the best way to include composition dependence yet.
      EquationOfStateOutputs<dim> eos_outputs (1);

      std::vector<double> viscosity_per_maxwell_arm(number_of_maxwell_arms, std::numeric_limits<double>::signaling_NaN());

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs);

          out.densities[i]                       = eos_outputs.densities[0];
          out.thermal_expansion_coefficients[i]  = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[i]                   = eos_outputs.specific_heat_capacities[0];
          out.compressibilities[i]               = eos_outputs.compressibilities[0];

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          for (unsigned int a=0; a<number_of_maxwell_arms; ++a)
            viscosity_per_maxwell_arm[a] =
              elastic_rheology.calculate_viscoelastic_viscosity(viscosities[a], elastic_shear_moduli[a]);

          // Parallel Maxwell arms would have the same strain rates so summing the stress contributions is equivalent to summing the viscosities.
          out.viscosities[i] = std::accumulate(viscosity_per_maxwell_arm.begin(), viscosity_per_maxwell_arm.end(), 0.0);
        }

      fill_elastic_outputs(in, elastic_shear_moduli, out);
      elastic_rheology.fill_elastic_additional_outputs(in, elastic_shear_moduli, out);
      elastic_rheology.fill_reaction_outputs(in, elastic_shear_moduli, out);
      elastic_rheology.fill_reaction_rates(in, elastic_shear_moduli, out);
    }



    template <int dim>
    double
    ViscoelasticComposite<dim>::elastic_timestep () const
    {
      double minimum_relaxation_time = std::numeric_limits<double>::max();

      for (unsigned int n = 0; n < number_of_maxwell_arms; ++n)
        {
          const double relaxation_time = viscosities[n] / elastic_shear_moduli[n];
          minimum_relaxation_time = std::min(minimum_relaxation_time, relaxation_time);
        }

      return minimum_relaxation_time;
    }




    template <int dim>
    void
    ViscoelasticComposite<dim>::fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                      const std::vector<double> &average_elastic_shear_moduli,
                                                      MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const std::shared_ptr<MaterialModel::ElasticOutputs<dim>> elastic_out = out.template get_additional_output_object<MaterialModel::ElasticOutputs<dim>>();
      const std::shared_ptr<HeatingModel::PrescribedShearHeatingOutputs<dim>>
      heating_out = out.template get_additional_output_object<HeatingModel::PrescribedShearHeatingOutputs<dim>>();

      if (!in.requests_property(MaterialProperties::additional_outputs))
        return;

      // The viscosity should be averaged if material averaging is applied.
      std::vector<double> effective_creep_viscosities;
      if (this->get_parameters().material_averaging != MaterialAveraging::none)
        {
          MaterialModelOutputs<dim> out_copy(out.n_evaluation_points(),
                                             this->introspection().n_compositional_fields);
          out_copy.viscosities = out.viscosities;

          const MaterialAveraging::AveragingOperation averaging_operation_for_viscosity = get_averaging_operation_for_viscosity(this->get_parameters().material_averaging);
          MaterialAveraging::average(averaging_operation_for_viscosity,
                                     in.current_cell,
                                     this->introspection().quadratures.velocities,
                                     this->get_mapping(),
                                     in.requested_properties,
                                     out_copy);

          effective_creep_viscosities = out_copy.viscosities;
        }
      else
        effective_creep_viscosities = out.viscosities;

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const SymmetricTensor<2,dim> deviatoric_strain_rate =
            Utilities::Tensors::consistent_deviator(in.strain_rate[i]);

          const double effective_creep_viscosity = effective_creep_viscosities[i];

          // Accumulate each arm's force contribution separately, then sum.
          std::vector<SymmetricTensor<2,dim>> elastic_force_per_arm(number_of_maxwell_arms);

          for (unsigned int n = 0; n < number_of_maxwell_arms; ++n)
            {
              const unsigned int stress_start_index_arm =
                this->introspection().compositional_index_for_name("ve_stress_xx_" + Utilities::int_to_string(n));

              const SymmetricTensor<2,dim> stress_0_advected (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index_arm],
                                                              &in.composition[i][stress_start_index_arm]+n_independent_components));

              const SymmetricTensor<2,dim> stress_old (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index_arm+n_independent_components],
                                                       &in.composition[i][stress_start_index_arm+n_independent_components]+n_independent_components));

              const double timestep_ratio = elastic_rheology.calculate_timestep_ratio();
              const double viscosity_ratio = elastic_rheology.calculate_viscoelastic_viscosity(viscosities[n], elastic_shear_moduli[n]) /
                                             elastic_rheology.calculate_elastic_viscosity(average_elastic_shear_moduli[n]);

              const SymmetricTensor<2,dim> stress = 2. * effective_creep_viscosity * deviatoric_strain_rate + viscosity_ratio * stress_0_advected
                                                    + (1. - timestep_ratio) * (1. - viscosity_ratio) * stress_old;

              const double dtc = timestep_ratio * this->elastic_timestep();

              elastic_force_per_arm[n] = -1. * (viscosity_ratio *stress_0_advected + (1. - timestep_ratio) * (1. - viscosity_ratio) * stress_old);

              // Assume incompressibility.
              const SymmetricTensor<2, dim> visco_plastic_strain_rate = deviatoric_strain_rate - ((stress - stress_0_advected) / (2. * dtc *average_elastic_shear_moduli[i]));
              // If compressible,
              // visco_plastic_strain_rate = visco_plastic_strain_rate -
              //                             1. / 3. * trace(visco_plastic_strain_rate) * unit_symmetric_tensor<dim>();

              // The shear heating term needs to account for the elastic stress, but only the visco_plastic strain rate.
              // This is best computed here, and stored for later use by the heating model.
              if (heating_out != nullptr)
                heating_out->prescribed_shear_heating_rates[i] = stress * visco_plastic_strain_rate;
            }

          if (elastic_out != nullptr)
            elastic_out->elastic_force[i] = std::accumulate(elastic_force_per_arm.begin(), elastic_force_per_arm.end(), SymmetricTensor<2,dim>());
        }
    }



    template <int dim>
    bool
    ViscoelasticComposite<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }



    template <int dim>
    void
    ViscoelasticComposite<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic composite");
        {
          prm.declare_entry ("Elastic shear moduli", "75.0e9",
                             Patterns::List(Patterns::Double(0.)),
                             "List of elastic shear moduli, $G$, for each of the maxwell arms. "
                             "The default value of 75 GPa is representative of mantle rocks. Units:\\si{\\pascal\\second}.");
          prm.declare_entry ("Viscosities", "1.e21",
                             Patterns::List(Patterns::Double (0.)),
                             "List of viscosities for each of the maxwell arms. "
                             "If only one value is given, then all use the same value. "
                             "Units: $\\text{Pa}\\text{s}$.");
          prm.declare_entry ("Number of maxwell arms", "2",
                             Patterns::Integer(0),
                             "Number of parallel maxwell arms in the composite model");

          prm.enter_subsection("Elastic parameters");
          Rheology::Elasticity<dim>::declare_parameters(prm);
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoelasticComposite<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Viscoelastic composite");
        {
          // For us the options mean different arms of the maxwell model, not different compositional fields
          const unsigned int number_of_maxwell_arms = prm.get_integer("Number of maxwell arms");


          viscosities = Utilities::possibly_extend_from_1_to_N (dealii::Utilities::string_to_double(Utilities::split_string_list(prm.get("Viscosities"))),
                                                                number_of_maxwell_arms, "Viscosities");
          elastic_shear_moduli = Utilities::possibly_extend_from_1_to_N (
                                   dealii::Utilities::string_to_double(
                                     Utilities::split_string_list(prm.get("Elastic shear moduli"))),
                                   number_of_maxwell_arms,
                                   "Elastic shear moduli");

          elastic_rheology.initialize_simulator (this->get_simulator());

          // Elastic rheology also has its own shear moduli, which we will
          // use when we have multiple compositional fields.
          prm.enter_subsection("Elastic parameters");
          elastic_rheology.parse_parameters(prm);
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ViscoelasticComposite<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      elastic_rheology.create_elastic_additional_outputs(out);
    }
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoelasticComposite,
                                   "viscoelastic composite",
                                   "A viscoelastic material model that uses a composite of parallel "
                                   "Maxwell arms to represent the viscoelastic behavior of the material. "
                                   "Each arm has its own viscosity and elastic shear modulus, allowing for "
                                   "a more complex representation of viscoelasticity.")
  }
}
