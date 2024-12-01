/*
  Copyright (C) 2019 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/rheology/elasticity.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace
    {
      std::vector<std::string> make_elastic_additional_outputs_names()
      {
        std::vector<std::string> names;
        names.emplace_back("elastic_shear_modulus");
        return names;
      }
    }

    template <int dim>
    ElasticAdditionalOutputs<dim>::ElasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_elastic_additional_outputs_names()),
      elastic_shear_moduli(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    std::vector<double>
    ElasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      (void)idx; // suppress warning in release mode
      AssertIndexRange (idx, 1);
      return elastic_shear_moduli;
    }



    namespace Rheology
    {
      template <int dim>
      void
      Elasticity<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Elastic shear moduli", "75.0e9",
                           Patterns::List(Patterns::Double(0.)),
                           "List of elastic shear moduli, $G$, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of all compositional fields or only "
                           "those corresponding to chemical compositions. "
                           "The default value of 75 GPa is representative of mantle rocks. Units: Pa.");
        prm.declare_entry ("Use fixed elastic time step", "unspecified",
                           Patterns::Selection("true|false|unspecified"),
                           "Select whether the material time scale in the viscoelastic constitutive "
                           "relationship uses the regular numerical time step or a separate fixed "
                           "elastic time step throughout the model run. The fixed elastic time step "
                           "is always used during the initial time step. If a fixed elastic time "
                           "step is used throughout the model run, a stress averaging scheme is "
                           "applied to account for differences with the numerical time step. An "
                           "alternative approach is to limit the maximum time step size so that it "
                           "is equal to the elastic time step. The default value of this parameter is "
                           "'unspecified', which throws an exception during runtime. In order for "
                           "the model to run the user must select 'true' or 'false'.");
        prm.declare_entry ("Fixed elastic time step", "1.e3",
                           Patterns::Double (0.),
                           "The fixed elastic time step $dte$. Units: years if the "
                           "'Use years in output instead of seconds' parameter is set; "
                           "seconds otherwise.");
        prm.declare_entry ("Stabilization time scale factor", "1.",
                           Patterns::Double (1.),
                           "A stabilization factor for the elastic stresses that influences how fast "
                           "elastic stresses adjust to deformation. This value is equal to the "
                           "elastic time step divided by the computational time step. "
                           "The default value of 1.0 may lead to oscillatory motion. "
                           "Increasing this factor to 2.0 can reduce oscillations while "
                           "preserving an immediate elastic response. In complex models the factor "
                           "can be increased further to improve convergence behaviour. "
                           "As the stabilization factor increases, the effective viscosity "
                           "gets smaller, and is balanced by an increasing body force term. "
                           "For composite rheologies that use this formulation of elasticity, "
                           "setting an infinite shear modulus only recovers the nonelastic part of "
                           "the rheology if this stabilization factor is equal to 1.0.");
        prm.declare_entry ("Elastic damper viscosity", "0.0",
                           Patterns::Double (0.),
                           "Viscosity of a viscous damper that acts in parallel with the elastic "
                           "element to stabilize behavior. Units: \\si{\\pascal\\second}");
      }



      template <int dim>
      void
      Elasticity<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Retrieve the list of composition names
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

        // Retrieve the list of names of fields that represent chemical compositions, and not, e.g.,
        // plastic strain
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();

        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(), "background");
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        Utilities::MapParsing::Options options(chemical_field_names, "Elastic shear moduli");
        options.list_of_allowed_keys = compositional_field_names;

        elastic_shear_moduli = Utilities::MapParsing::parse_map_to_double_array (prm.get("Elastic shear moduli"),
                                                                                 options);

        // Stabilize elasticity through a viscous damper
        elastic_damper_viscosity = prm.get_double("Elastic damper viscosity");

        if (prm.get ("Use fixed elastic time step") == "true")
          use_fixed_elastic_time_step = true;
        else if (prm.get ("Use fixed elastic time step") == "false")
          use_fixed_elastic_time_step = false;
        else
          AssertThrow(false, ExcMessage("'Use fixed elastic time step' must be set to 'true' or 'false'"));

        stabilization_time_scale_factor = prm.get_double ("Stabilization time scale factor");

        fixed_elastic_time_step = prm.get_double ("Fixed elastic time step");
        AssertThrow(fixed_elastic_time_step > 0,
                    ExcMessage("The fixed elastic time step must be greater than zero"));

        if (this->convert_output_to_years())
          fixed_elastic_time_step *= year_in_seconds;

        AssertThrow(this->get_parameters().enable_elasticity == true,
                    ExcMessage ("Material model Viscoelastic only works if 'Enable elasticity' is set to true"));

        // Check whether the compositional fields representing the viscoelastic
        // stress tensor are both named correctly and listed in the right order.
        AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xx") == 0,
                    ExcMessage("Rheology model Elasticity only works if the first "
                               "compositional field is called ve_stress_xx."));
        AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yy") == 1,
                    ExcMessage("Rheology model Elasticity only works if the second "
                               "compositional field is called ve_stress_yy."));
        if (dim == 2)
          {
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xy") == 2,
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field is called ve_stress_xy."));
          }
        else if (dim == 3)
          {
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_zz") == 2,
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field is called ve_stress_zz."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xy") == 3,
                        ExcMessage("Rheology model Elasticity only works if the fourth "
                                   "compositional field is called ve_stress_xy."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_xz") == 4,
                        ExcMessage("Rheology model Elasticity only works if the fifth "
                                   "compositional field is called ve_stress_xz."));
            AssertThrow(this->introspection().compositional_index_for_name("ve_stress_yz") == 5,
                        ExcMessage("Rheology model Elasticity only works if the sixth "
                                   "compositional field is called ve_stress_yz."));
          }
        else
          AssertThrow(false, ExcNotImplemented());

        // Currently, it only makes sense to use this material model when the nonlinear solver
        // scheme does a single Advection iteration and at minimum one Stokes iteration. More
        // than one nonlinear Advection iteration will produce an unrealistic build-up of
        // viscoelastic stress, which are tracked through compositional fields.
        AssertThrow((this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                     ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::single_Advection_iterated_Stokes
                     ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes
                     ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::single_Advection_iterated_defect_correction_Stokes),
                    ExcMessage("The material model will only work with the nonlinear "
                               "solver schemes 'single Advection, single Stokes', "
                               "'single Advection, iterated Stokes', "
                               "'single Advection, iterated Newton Stokes', and "
                               "'single Advection, iterated defect correction Stokes' "));

        // Functionality to average the additional RHS terms over the cell is not implemented.
        // Also, there is no option implemented in this rheology module to project to Q1 the viscosity
        // in the elastic force term for the RHS.
        // Consequently, it is only possible to use elasticity with the Material averaging schemes
        // 'none', 'harmonic average only viscosity', and 'geometric average only viscosity'.
        // TODO: Find a way to include 'project to Q1 only viscosity'.
        AssertThrow((this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::none
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::harmonic_average_only_viscosity
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::geometric_average_only_viscosity
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::default_averaging),
                    ExcMessage("Material models with elasticity can only be used with the material "
                               "averaging schemes 'none', 'harmonic average only viscosity' and "
                               "'geometric average only viscosity'. This parameter ('Material averaging') "
                               "is located within the 'Material model' subsection."));
      }



      template <int dim>
      void
      Elasticity<dim>::create_elastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<ElasticAdditionalOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<ElasticAdditionalOutputs<dim>> (n_points));
          }
      }



      template <int dim>
      void
      Elasticity<dim>::fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                             const std::vector<double> &average_elastic_shear_moduli,
                                             MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic outputs
        MaterialModel::ElasticOutputs<dim>
        *elastic_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();

        if (elastic_out == nullptr)
          return;

        if (in.requests_property(MaterialProperties::additional_outputs))
          {
            // The viscosity should be averaged if material averaging is applied.
            std::vector<double> effective_creep_viscosities;
            if (this->get_parameters().material_averaging != MaterialAveraging::none)
              {
                MaterialModelOutputs<dim> out_copy(out.n_evaluation_points(),
                                                   this->introspection().n_compositional_fields);
                out_copy.viscosities = out.viscosities;

                const MaterialAveraging::AveragingOperation averaging_operation_for_viscosity =
                  get_averaging_operation_for_viscosity(this->get_parameters().material_averaging);
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
                // Get old stresses from compositional fields
                const SymmetricTensor<2,dim> stress_old (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][0],
                                                         &in.composition[i][0]+SymmetricTensor<2,dim>::n_independent_components));

                elastic_out->elastic_force[i] = -effective_creep_viscosities[i] / calculate_elastic_viscosity(average_elastic_shear_moduli[i]) * stress_old;
                // The viscoelastic strain rate is needed only when the Newton method is selected.
                const typename Parameters<dim>::NonlinearSolver::Kind nonlinear_solver = this->get_parameters().nonlinear_solver;
                if ((nonlinear_solver == Parameters<dim>::NonlinearSolver::iterated_Advection_and_Newton_Stokes) ||
                    (nonlinear_solver == Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes))
                  elastic_out->viscoelastic_strain_rate[i] = calculate_viscoelastic_strain_rate(in.strain_rate[i], stress_old, average_elastic_shear_moduli[i]);
              }
          }
      }



      template <int dim>
      void
      Elasticity<dim>::fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                              const std::vector<double> &average_elastic_shear_moduli,
                                              MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
          {
            // Get old (previous time step) velocity gradients
            std::vector<Point<dim>> quadrature_positions(in.n_evaluation_points());
            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

            std::vector<double> solution_values(this->get_fe().dofs_per_cell);
            in.current_cell->get_dof_values(this->get_old_solution(),
                                            solution_values.begin(),
                                            solution_values.end());

            // Only create the evaluator the first time we get here
            if (!evaluator)
              evaluator = std::make_unique<FEPointEvaluation<dim,dim>>(this->get_mapping(),
                                                                        this->get_fe(),
                                                                        update_gradients,
                                                                        this->introspection().component_indices.velocities[0]);

            // Initialize the evaluator for the old velocity gradients
            evaluator->reinit(in.current_cell, quadrature_positions);
            evaluator->evaluate(solution_values,
                                EvaluationFlags::gradients);

            const double dte = elastic_timestep();
            const double dt = this->get_timestep();

            // The viscosity should be averaged if material averaging is applied.
            // Here the averaging scheme "project to Q1 (only viscosity)"  is
            // excluded, because there is no way to know the quadrature formula
            // used for evaluation.
            // TODO: find a way to include "project to Q1 (only viscosity)" as well.
            std::vector<double> effective_creep_viscosities;
            if (this->get_parameters().material_averaging != MaterialAveraging::none &&
                this->get_parameters().material_averaging != MaterialAveraging::project_to_Q1 &&
                this->get_parameters().material_averaging != MaterialAveraging::project_to_Q1_only_viscosity)
              {
                MaterialModelOutputs<dim> out_copy(out.n_evaluation_points(),
                                                   this->introspection().n_compositional_fields);
                out_copy.viscosities = out.viscosities;

                const MaterialAveraging::AveragingOperation averaging_operation_for_viscosity =
                  get_averaging_operation_for_viscosity(this->get_parameters().material_averaging);
                MaterialAveraging::average(averaging_operation_for_viscosity,
                                           in.current_cell,
                                           Quadrature<dim>(quadrature_positions),
                                           this->get_mapping(),
                                           in.requested_properties,
                                           out_copy);

                effective_creep_viscosities = out_copy.viscosities;
              }
            else
              effective_creep_viscosities = out.viscosities;

            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                // Get old stresses from compositional fields
                const SymmetricTensor<2,dim> stress_old(Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][0],
                                                        &in.composition[i][0]+SymmetricTensor<2,dim>::n_independent_components));

                // Calculate the rotated stresses
                // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
                const Tensor<2,dim> rotation = 0.5 * (evaluator->get_gradient(i) - transpose(evaluator->get_gradient(i)));

                // Calculate the current (new) stored elastic stress, which is a function of the material
                // properties (viscoelastic viscosity, shear modulus), elastic time step size, strain rate,
                // vorticity, prior (inherited) viscoelastic stresses and viscosity of the elastic damper.
                // In the absence of the elastic damper, the expression for "stress_new" is identical
                // to the one found in Moresi et al. (2003, J. Comp. Phys., equation 29).
                const double damped_elastic_viscosity = calculate_elastic_viscosity(average_elastic_shear_moduli[i]);

                // stress_0 is the combination of the elastic stress tensor stored at the end of the last time step and the change in that stress generated by local rotation
                const SymmetricTensor<2,dim> stress_0 = (stress_old + dte * ( symmetrize(rotation * Tensor<2,dim>(stress_old) ) - symmetrize(Tensor<2,dim>(stress_old) * rotation) ) );

                // stress_creep is the stress experienced by the viscous and elastic components.
                Assert(std::isfinite(in.strain_rate[i].norm()),
                       ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                  "not filled by the caller."));
                const SymmetricTensor<2,dim> stress_creep = 2. * effective_creep_viscosities[i] * ( deviator(in.strain_rate[i]) + stress_0 / (2. * damped_elastic_viscosity ) );

                // stress_new is the (new) stored elastic stress
                SymmetricTensor<2,dim> stress_new = stress_creep * (1. - (elastic_damper_viscosity / damped_elastic_viscosity)) + elastic_damper_viscosity * stress_0 / damped_elastic_viscosity;

                // Stress averaging scheme to account for difference between the elastic time step
                // and the numerical time step (see equation 32 in Moresi et al., 2003, J. Comp. Phys.)
                // Note that if there is no difference between the elastic timestep and the numerical
                // timestep, then no averaging occurs as dt/dte = 1.
                stress_new = ( ( 1. - ( dt / dte ) ) * stress_old ) + ( ( dt / dte ) * stress_new ) ;

                const SymmetricTensor<2,dim> stress_update = stress_new - stress_old;

                // Fill reaction terms
                Utilities::Tensors::unroll_symmetric_tensor_into_array(stress_update,
                                                                       &out.reaction_terms[i][0],
                                                                       &out.reaction_terms[i][0]+SymmetricTensor<2,dim>::n_independent_components);
              }
          }
      }



      template <int dim>
      double
      Elasticity<dim>::elastic_timestep () const
      {
        // The elastic time step (dte) is equal to the numerical time step if the time step number
        // is greater than 0 and the parameter 'use_fixed_elastic_time_step' is set to false.
        // On the first (0) time step the elastic time step is always equal to the value
        // specified in 'fixed_elastic_time_step', which is also used in all subsequent time
        // steps if 'use_fixed_elastic_time_step' is set to true.
        //
        // We also use this parameter when we are still *before* the first time step,
        // i.e., if the time step number is numbers::invalid_unsigned_int.
        const double dte = ( ( this->get_timestep_number() > 0 &&
                               this->simulator_is_past_initialization() &&
                               use_fixed_elastic_time_step == false )
                             ?
                             this->get_timestep() * stabilization_time_scale_factor
                             :
                             fixed_elastic_time_step);
        return dte;
      }



      template <int dim>
      const std::vector<double> &
      Elasticity<dim>::get_elastic_shear_moduli () const
      {
        return elastic_shear_moduli;
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_elastic_viscosity (const double shear_modulus) const
      {
        return shear_modulus*elastic_timestep() + elastic_damper_viscosity;
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_viscoelastic_viscosity (const double viscosity,
                                        const double shear_modulus) const
      {
        const double elastic_viscosity = calculate_elastic_viscosity(shear_modulus);
        return 1. / (1./elastic_viscosity + 1./viscosity);
      }



      template <int dim>
      SymmetricTensor<2,dim>
      Elasticity<dim>::
      calculate_viscoelastic_strain_rate(const SymmetricTensor<2,dim> &strain_rate,
                                         const SymmetricTensor<2,dim> &stored_stress,
                                         const double shear_modulus) const
      {
        // The first term in the following expression is the deviator of the true strain rate
        // of one or more isostress rheological elements (in series).
        // One of these elements must be an elastic component (potentially damped).
        // The second term corresponds to a fictional strain rate arising from
        // elastic stresses stored from the last time step.
        // Note the parallels with the viscous part of the strain rate deviator,
        // which is equal to 0.5 * stress / viscosity.
        return deviator(strain_rate) + 0.5 * deviator(stored_stress) /
               calculate_elastic_viscosity(shear_modulus);
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class ElasticAdditionalOutputs<dim>; \
  \
  namespace Rheology \
  { \
    template class Elasticity<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
