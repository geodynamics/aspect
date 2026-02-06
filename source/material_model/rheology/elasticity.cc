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

#include <aspect/utilities.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/material_model/viscoelastic.h>
#include <aspect/heating_model/shear_heating.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
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
        // We do not output all the additional output to
        // visualization output, only the ones listed here.
        names.emplace_back("elastic_shear_modulus");
        names.emplace_back("elastic_viscosity");
        return names;
      }
    }

    template <int dim>
    ElasticAdditionalOutputs<dim>::ElasticAdditionalOutputs (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_elastic_additional_outputs_names()),
      elastic_shear_moduli(n_points, numbers::signaling_nan<double>()),
      elastic_viscosity(n_points, numbers::signaling_nan<double>()),
      deviatoric_stress(n_points, SymmetricTensor<2,dim>())
    {}



    template <int dim>
    std::vector<double>
    ElasticAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
    {
      (void)idx; // suppress warning in release mode
      AssertIndexRange (idx, 2);

      switch (idx)
        {
          case 0:
            return elastic_shear_moduli;

          case 1:
            return elastic_viscosity;

          default:
            AssertThrow(false, ExcInternalError());
        }

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
                           "The fixed elastic time step $dte$. It is always used during the first "
                           "timestep; afterwards on if 'Used fixed elastic time step' is true. "
                           "Units: years if the 'Use years instead of seconds' parameter is set; "
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
                           "element to stabilize behavior. Units: $\\text{Pa}\\text{s}$");
      }



      template <int dim>
      void
      Elasticity<dim>::parse_parameters (ParameterHandler &prm)
      {
        AssertThrow(this->get_parameters().enable_elasticity == true,
                    ExcMessage("Rheology model elasticity only works if 'Enable elasticity' is set to true"));

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

        // When using the visco_plastic or viscoelastic material model,
        // make sure that no damping is applied. Damping could potentially
        // improve stability under rapidly changing dynamics, but
        // so far it has not been necessary.
        AssertThrow(elastic_damper_viscosity == 0.,
                    ExcMessage("The viscoelastic material model and the visco-plastic material model with elasticity enabled require "
                               "that no elastic damping is applied."));
        // An update of the stored stresses is done in an operator splitting step for fields or by the particle property 'elastic stress'.
        AssertThrow(this->get_parameters().use_operator_splitting || (this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx")),
                    ExcMessage("The viscoelastic material model and the visco-plastic material model with elasticity enabled require "
                               "operator splitting for stresses tracked on compositional fields or the particle property 'elastic stress' "
                               "for stresses tracked on particles."));
        // If the operator splitting scheme is used, make sure to use its fixed step solver, as we know the update and it should be applied in one step.
        if (this->get_parameters().use_operator_splitting)
          AssertThrow(this->get_parameters().reaction_solver_type == Parameters<dim>::ReactionSolverType::fixed_step,
                      ExcMessage("If the operator splitting scheme is used, its solver should be set to 'fixed step'."));
        if ((this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx")))
          AssertThrow(!this->get_parameters().use_operator_splitting,
                      ExcMessage("If stresses are tracked on particles, the stress update is applied by the particle property 'elastic stress' "
                                 "and operator splitting should not be turned on. "));

        // Check that 3+3 in 2D or 6+6 in 3D stress fields exist.
        AssertThrow((this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::stress) == 2*SymmetricTensor<2,dim>::n_independent_components),
                    ExcMessage("Rheology model Elasticity requires 3+3 in 2D or 6+6 in 3D fields of type stress."));

        // Check that the compositional fields representing the viscoelastic
        // stress tensor components are both named correctly and listed in the right order
        // as well as use a discontinuous discretization.
        std::vector<std::string> stress_field_names = this->introspection().get_names_for_fields_of_type(CompositionalFieldDescription::stress);
        std::vector<unsigned int> stress_field_indices = this->introspection().get_indices_for_fields_of_type(CompositionalFieldDescription::stress);

        // The discontinuous element is required to accommodate discontinuous
        // strain rates that feed into the stored stresses.
        const std::vector<bool> use_discontinuous_composition_discretization = this->get_parameters().use_discontinuous_composition_discretization;
        for (auto stress_field : stress_field_indices)
          AssertThrow(use_discontinuous_composition_discretization[stress_field],
                      ExcMessage("The viscoelastic material model and the visco-plastic material model with elasticity enabled require "
                                 "the use of discontinuous elements for compositions that represent stress tensor components."));

        // We require a consecutive range of indices (for example for FEPointEvaluation)
        // to extract the fields representing the viscoelastic stress tensor components,
        // so check that they are listed without interruption by other fields.
        // They do not, however, have to be the first fields listed.
        AssertThrow(((stress_field_indices[2*n_independent_components-1] - stress_field_indices[0]) == (2*n_independent_components-1)),
                    ExcMessage("Rheology model Elasticity requires that the compositional fields representing stress tensor components are listed in consecutive order."));

        AssertThrow(stress_field_names[0] == "ve_stress_xx",
                    ExcMessage("Rheology model Elasticity only works if the first "
                               "compositional field representing stress tensor components is called ve_stress_xx."));
        AssertThrow(stress_field_names[1] == "ve_stress_yy",
                    ExcMessage("Rheology model Elasticity only works if the second "
                               "compositional field representing stress tensor components is called ve_stress_yy."));
        if (dim == 2)
          {
            AssertThrow(stress_field_names[2] == "ve_stress_xy",
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field representing stress tensor components is called ve_stress_xy."));

            AssertThrow(stress_field_names[3] == "ve_stress_xx_old",
                        ExcMessage("Rheology model Elasticity only works if the fourth "
                                   "compositional field representing stress tensor components is called ve_stress_xx_old."));
            AssertThrow(stress_field_names[4] == "ve_stress_yy_old",
                        ExcMessage("Rheology model Elasticity only works if the fifth "
                                   "compositional field representing stress tensor components is called ve_stress_yy_old."));
            AssertThrow(stress_field_names[5] == "ve_stress_xy_old",
                        ExcMessage("Rheology model Elasticity only works if the sixth "
                                   "compositional field representing stress tensor components is called ve_stress_xy_old."));
          }
        else if (dim == 3)
          {
            AssertThrow(stress_field_names[2] == "ve_stress_zz",
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field representing stress tensor components is called ve_stress_zz."));
            AssertThrow(stress_field_names[3] == "ve_stress_xy",
                        ExcMessage("Rheology model Elasticity only works if the fourth "
                                   "compositional field representing stress tensor components is called ve_stress_xy."));
            AssertThrow(stress_field_names[4] == "ve_stress_xz",
                        ExcMessage("Rheology model Elasticity only works if the fifth "
                                   "compositional field representing stress tensor components is called ve_stress_xz."));
            AssertThrow(stress_field_names[5] == "ve_stress_yz",
                        ExcMessage("Rheology model Elasticity only works if the sixth "
                                   "compositional field representing stress tensor components is called ve_stress_yz."));

            AssertThrow(stress_field_names[6] == "ve_stress_xx_old",
                        ExcMessage("Rheology model Elasticity only works if the seventh "
                                   "compositional field representing stress tensor components is called ve_stress_xx_old."));
            AssertThrow(stress_field_names[7] == "ve_stress_yy_old",
                        ExcMessage("Rheology model Elasticity only works if the eighth "
                                   "compositional field representing stress tensor components is called ve_stress_yy_old."));
            AssertThrow(stress_field_names[8] == "ve_stress_zz_old",
                        ExcMessage("Rheology model Elasticity only works if the ninth "
                                   "compositional field representing stress tensor components is called ve_stress_zz_old."));
            AssertThrow(stress_field_names[9] == "ve_stress_xy_old",
                        ExcMessage("Rheology model Elasticity only works if the tenth "
                                   "compositional field representing stress tensor components is called ve_stress_xy_old."));
            AssertThrow(stress_field_names[10] == "ve_stress_xz_old",
                        ExcMessage("Rheology model Elasticity only works if the eleventh "
                                   "compositional field representing stress tensor components is called ve_stress_xz_old."));
            AssertThrow(stress_field_names[11] == "ve_stress_yz_old",
                        ExcMessage("Rheology model Elasticity only works if the twelfth "
                                   "compositional field representing stress tensor components is called ve_stress_yz_old."));
          }
        else
          AssertThrow(false, ExcNotImplemented());

        // We need to iterate over the Advection and Stokes equations.
        AssertThrow((this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::iterated_Advection_and_Stokes ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
                     this->get_parameters().nonlinear_solver ==
                     Parameters<dim>::NonlinearSolver::iterated_Advection_and_defect_correction_Stokes),
                    ExcMessage("The material model will only work with the nonlinear "
                               "solver schemes 'iterated Advection and Stokes', "
                               "'iterated Advection and defect correction Stokes', "
                               "'iterated Advection and Newton Stokes'."));

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
      Elasticity<dim>::create_elastic_additional_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create the ElasticAdditionalOutputs that include the average shear modulus, elastic
        // viscosity, timestep ratio and total deviatoric stress of the current timestep.
        if (out.template has_additional_output_object<ElasticAdditionalOutputs<dim>>() == false)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<ElasticAdditionalOutputs<dim>> (n_points));
          }

        // We need to modify the shear heating outputs to correctly account for elastic stresses.
        if (out.template has_additional_output_object<HeatingModel::PrescribedShearHeatingOutputs<dim>>() == false)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<HeatingModel::PrescribedShearHeatingOutputs<dim>> (n_points));
          }

        // Create the ReactionRateOutputs that are necessary for the operator splitting
        // step (either on the fields or directly on the particles)
        // that sets both sets of stresses to the total stress of the
        // previous timestep.
        if (out.template has_additional_output_object<ReactionRateOutputs<dim>>() == false &&
            (this->get_parameters().use_operator_splitting || (this->get_parameters().mapped_particle_properties).count(this->introspection().compositional_index_for_name("ve_stress_xx"))))
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std::make_unique<MaterialModel::ReactionRateOutputs<dim>>(n_points, this->n_compositional_fields()));
          }
      }



      template <int dim>
      void
      Elasticity<dim>::fill_elastic_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                             const std::vector<double> &average_elastic_shear_moduli,
                                             MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic outputs.
        // The structure is created during the Stokes assembly.
        const std::shared_ptr<MaterialModel::ElasticOutputs<dim>>
        elastic_out = out.template get_additional_output_object<MaterialModel::ElasticOutputs<dim>>();

        // Create a reference to the structure for the prescribed shear heating outputs.
        // The structure is created during the advection assembly.
        const std::shared_ptr<HeatingModel::PrescribedShearHeatingOutputs<dim>>
        heating_out = out.template get_additional_output_object<HeatingModel::PrescribedShearHeatingOutputs<dim>>();

        if (elastic_out == nullptr && heating_out == nullptr)
          return;

        // TODO should a RHS term be a separate MaterialProperties?
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

            const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");

            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                const SymmetricTensor<2, dim> deviatoric_strain_rate = Utilities::Tensors::consistent_deviator(in.strain_rate[i]);

                // Get stress from timestep $t$ rotated and advected into the current
                // timestep $t+\Delta t_c$ from the compositional fields.
                // This function is only evaluated during the assembly of the Stokes equations
                // (the force term goes into the rhs of the momentum equation).
                // This happens after the advection equations have been solved, and hence in.composition
                // contains the rotated and advected stresses $tau^{0adv}$.
                // Only at the beginning of the next timestep do we add the stress update of the
                // current timestep to the stress stored in the compositional fields, giving
                // $\tau{t+\Delta t_c}$ with $t+\Delta t_c$ being the current timestep.
                const SymmetricTensor<2,dim> stress_0_advected (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                                &in.composition[i][stress_start_index]+n_independent_components));

                // Get the old stress that is used to interpolate to timestep $t+\Delta t_c$. It is stored on the
                // second set of n_independent_components fields, e.g. in 2D on field 3, 4 and 5.
                const SymmetricTensor<2,dim> stress_old (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index+n_independent_components],
                                                         &in.composition[i][stress_start_index+n_independent_components]+n_independent_components));

                // Average effective creep viscosity
                // Use the viscosity corresponding to the stresses selected above.
                // out.viscosities is computed during the assembly of the Stokes equations
                // based on the current_linearization_point. This means that it will be updated after every
                // nonlinear Stokes iteration.
                // The effective creep viscosity has already been scaled with the timestep ratio dtc/dte.
                const double effective_creep_viscosity = effective_creep_viscosities[i];

                // The force term is computed as:
                // $\frac{-\eta_{effcreep} \tau_{0adv}}{\eta_{e}}$, where $\eta_{effcreep}$ is the
                // current harmonic average of the viscous and elastic viscosity, or the yield stress
                // divided by two times the second invariant of the deviatoric strain rate.
                // In case the computational timestep differs from the elastic timestep,
                // linearly interpolate between the two.
                const double timestep_ratio = calculate_timestep_ratio();
                // The elastic viscosity has also already been scaled with the timestep ratio.
                const double viscosity_ratio = effective_creep_viscosity / calculate_elastic_viscosity(average_elastic_shear_moduli[i]);

                if (elastic_out != nullptr)
                  {
                    elastic_out->elastic_force[i] = -1. * (viscosity_ratio * stress_0_advected
                                                           + (1. - timestep_ratio) * (1. - viscosity_ratio) * stress_old);

                    // The viscoelastic strain rate is needed only when the Newton method is selected.
                    const typename Parameters<dim>::NonlinearSolver::Kind nonlinear_solver = this->get_parameters().nonlinear_solver;
                    if ((nonlinear_solver == Parameters<dim>::NonlinearSolver::iterated_Advection_and_Newton_Stokes) ||
                        (nonlinear_solver == Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes))
                      elastic_out->viscoelastic_strain_rate[i] = calculate_viscoelastic_strain_rate(
                                                                   in.strain_rate[i], stress_0_advected, stress_old, effective_creep_viscosity, average_elastic_shear_moduli[i]);
                  }

                // Apply the stress update to get the total stress of timestep t.
                const SymmetricTensor<2, dim> stress = 2. * effective_creep_viscosity * deviatoric_strain_rate + viscosity_ratio * stress_0_advected +
                                                       (1. - timestep_ratio) * (1. - viscosity_ratio) * stress_old;

                // Obtain the computational timestep by multiplying the ratio between the computational
                // and elastic timestep $\frac{\Delta t_c}{\Delta t_{el}}$ with the elastic timestep.
                const double dtc = timestep_ratio * elastic_timestep();

                // Assume incompressibility.
                const SymmetricTensor<2, dim> visco_plastic_strain_rate = deviatoric_strain_rate - ((stress - stress_0_advected) / (2. * dtc * average_elastic_shear_moduli[i]));
                // If compressible,
                // visco_plastic_strain_rate = visco_plastic_strain_rate -
                //                             1. / 3. * trace(visco_plastic_strain_rate) * unit_symmetric_tensor<dim>();

                // The shear heating term needs to account for the elastic stress, but only the visco_plastic strain rate.
                // This is best computed here, and stored for later use by the heating model.
                if (heating_out != nullptr)
                  heating_out->prescribed_shear_heating_rates[i] = stress * visco_plastic_strain_rate;
              }

          }
      }




      template <int dim>
      void
      Elasticity<dim>::fill_elastic_additional_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                        const std::vector<double> &average_elastic_shear_moduli,
                                                        MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic additional outputs
        const std::shared_ptr<MaterialModel::ElasticAdditionalOutputs<dim>>
        elastic_additional_out = out.template get_additional_output_object<MaterialModel::ElasticAdditionalOutputs<dim>>();

        if (elastic_additional_out == nullptr || !in.requests_property(MaterialProperties::additional_outputs))
          return;

        const double timestep_ratio = calculate_timestep_ratio();
        const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");

        for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
          {
            const double eta = out.viscosities[i];
            const SymmetricTensor<2, dim> deviatoric_strain_rate = Utilities::Tensors::consistent_deviator(in.strain_rate[i]);
            const SymmetricTensor<2,dim> stress_0_advected (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                            &in.composition[i][stress_start_index]+n_independent_components));
            const SymmetricTensor<2,dim> stress_old (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index+n_independent_components],
                                                     &in.composition[i][stress_start_index+n_independent_components]+n_independent_components));

            const double elastic_viscosity = calculate_elastic_viscosity(average_elastic_shear_moduli[i]);

            // Apply the stress update to get the total deviatoric stress of timestep t.
            elastic_additional_out->deviatoric_stress[i] = 2. * eta * deviatoric_strain_rate + eta / elastic_viscosity * stress_0_advected + (1. - timestep_ratio) * (1. - eta / elastic_viscosity) * stress_old;
            elastic_additional_out->elastic_shear_moduli[i] = average_elastic_shear_moduli[i];
            elastic_additional_out->elastic_viscosity[i] = elastic_viscosity;
          }
      }



      // Rotate the stresses of the previous timestep $t$ into the current timestep $t+dtc$.
      template <int dim>
      void
      Elasticity<dim>::fill_reaction_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                              const std::vector<double> &,
                                              MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (in.current_cell.state() == IteratorState::valid
            && in.requests_property(MaterialProperties::reaction_terms))
          {
            // Get the velocity gradients of the current timestep $t+dtc$
            // at the requested location in in.position.
            std::vector<Point<dim>> quadrature_positions(in.n_evaluation_points());
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

            // Get the current velocity gradients, which get
            // updated in each nonlinear iteration.
            // This means we use the rotation tensor W^(t+dtc), not W^(t).
            std::vector<double> solution_values(this->get_fe().dofs_per_cell);
            in.current_cell->get_dof_values(this->get_current_linearization_point(),
                                            solution_values.begin(),
                                            solution_values.end());

            // Only create the evaluator the first time we get here.
            if (!evaluator)
              evaluator = std::make_unique<FEPointEvaluation<dim,dim>>(this->get_mapping(),
                                                                        this->get_fe(),
                                                                        update_gradients,
                                                                        this->introspection().component_indices.velocities[0]);

            // Initialize the evaluator for the velocity gradients.
            evaluator->reinit(in.current_cell, quadrature_positions);
            evaluator->evaluate(solution_values,
                                EvaluationFlags::gradients);

            // Get the fully updated stress of the previous timestep $t$.
            const std::vector<SymmetricTensor<2, dim>> stress_t = retrieve_stress_previous_timestep(in, quadrature_positions);

            Assert(out.reaction_terms.size() == in.n_evaluation_points(), ExcMessage("Out reaction terms not equal to n eval points."));

            const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");

            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                Assert(out.reaction_terms[i].size() == this->n_compositional_fields(), ExcMessage("Out reaction terms i not equal to n fields."));

                // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
                const Tensor<2, dim> rotation = 0.5 * (evaluator->get_gradient(i) - transpose(evaluator->get_gradient(i)));

                // stress_0 (i.e., $\tau^{0}$) is the sum of the stress tensor stored at the end of the last time step (stress_t)
                // and the change in that stress generated by local rotation over the computational timestep $\Delta t_c$:
                // stress_0 = stress_t + dtc * (symmetrize(rotation * stress_t) - symmetrize(stress_t * rotation)).
                // The reaction terms should be filled with the change from the previous stress to the rotated stress_0.
                // In case fields are used, in.composition holds the current_linearization_point,
                // which for the first nonlinear iteration means the extrapolated solution from
                // the last and previous to last timesteps. In later iterations, it holds the current solution.
                // If we subtract the stress in in.composition from stress_0 to compute the stress change, we would alternately compute
                // a reaction term that is correct and one that is zero. Therefore we subtract the old stress stress_t we retrieved
                // above. This means the reaction terms are:
                // stress_0 - stress_t = stress_t + dtc * (symmetrize(rotation * stress_t) - symmetrize(stress_t * rotation)) - stress_t
                // = dtc * (symmetrize(rotation * stress_t) - symmetrize(stress_t * rotation)).
                // TODO make sure these are equivalent and then rm
                //const SymmetricTensor<2, dim> stress_change = this->get_timestep() * (symmetrize(rotation * Tensor<2, dim>(stress_t[i])) - symmetrize(Tensor<2, dim>(stress_t[i]) * rotation));
                const SymmetricTensor<2, dim> stress_change = this->get_timestep() * symmetrize(rotation * Tensor<2, dim>(stress_t[i]) - Tensor<2, dim>(stress_t[i]) * rotation);

                Utilities::Tensors::unroll_symmetric_tensor_into_array(stress_change,
                                                                       &out.reaction_terms[i][stress_start_index],
                                                                       &out.reaction_terms[i][stress_start_index]+n_independent_components);
              }
          }
      }

      // The following function computes the reaction rates for the operator
      // splitting step that at the beginning of the new timestep $t+dtc$ updates the
      // stored compositions $tau^{0\mathrm{adv}}$ at time $t$ to $tau^{t}$.
      // This update consists of the stress change resulting from system evolution,
      // but does not advect or rotate the stress tensor. Advection is done by
      // solving the advection equation and the stress tensor is rotated through
      // the source term (reaction_terms) of that same equation.
      template <int dim>
      void
      Elasticity<dim>::fill_reaction_rates (const MaterialModel::MaterialModelInputs<dim> &in,
                                            const std::vector<double> &average_elastic_shear_moduli,
                                            MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        const std::shared_ptr<ReactionRateOutputs<dim>> reaction_rate_out
          = out.template get_additional_output_object<ReactionRateOutputs<dim>>();

        if (reaction_rate_out == nullptr)
          return;

        // Set all reaction rates to zero
        // TODO Should this only set those rates to zero
        // that are used to update the stresses instead of all rates?
        // What if other rheologies also fill reaction rates?
        for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
          for (unsigned int c = 0; c < in.composition[i].size(); ++c)
            reaction_rate_out->reaction_rates[i][c] = 0.0;

        // It doesn't make sense to update the stresses at time zero;
        // return reaction rates of zero.
        if (this->get_timestep_number() == 0)
          return;

        // At the moment when the reaction rates are required (at the beginning of the timestep),
        // the solution vector 'solution' holds the stress from the previous timestep,
        // advected into the new position of the previous timestep, so $\tau^{t}_{0adv}$.
        // This is the same as the vector 'old_solution' holds. At later moments during the current timestep,
        // 'solution' will hold the current_linearization_point instead of the solution of the previous timestep.
        //
        // In case fields are used to track the stresses, MaterialModelInputs are based on 'solution'
        // when calling the MaterialModel for the reaction rates. When particles are used, MaterialModelInputs
        // for this function are filled with the old_solution (including for the strain rate), except for the
        // compositions that represent the stress tensor components, these are taken directly from the
        // particles in the property plugin by default (although this can be changed from the input file).
        // As the particles are restored to their pre-advection location at the beginning of each nonlinear iteration,
        // their values and positions correspond to the old solution.
        // This means that in both cases we can use 'in' to get to the $\tau^{t}_{0adv}$ and velocity/strain rate of the
        // previous timestep.
        // TODO The additional outputs include the reaction rates, so we also have to fill the reaction_rates
        // if additional_outputs are requested.
        if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 &&
            (in.requests_property(MaterialProperties::reaction_rates) || in.requests_property(MaterialProperties::additional_outputs)))
          {
            const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");

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
                                           this->introspection().quadratures.velocities,
                                           this->get_mapping(),
                                           in.requested_properties,
                                           out_copy);

                effective_creep_viscosities = out_copy.viscosities;
              }
            else
              effective_creep_viscosities = out.viscosities;


            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                // Get $\tau^{0adv}$ of the previous timestep t from the compositional fields.
                // This stress includes the rotation and advection of the previous timestep,
                // i.e., the reaction term (which prescribes the change in stress due to rotation
                // over the previous timestep) has already been applied during the previous timestep.
                const SymmetricTensor<2, dim> stress_0_t (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                          &in.composition[i][stress_start_index]+n_independent_components));

                // Get the old stress that is used to interpolate to timestep $t+\Delta t_c$. It is stored on the
                // second set of n_independent_components fields, e.g. in 2D on field 3, 4 and 5.
                // The old stress was advected into the previous timestep, but not rotated.
                // Below we update it to full stress of the previous timestep, so that it can be
                // advected into the current timestep.
                const SymmetricTensor<2, dim>  stress_old (Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index+n_independent_components],
                                                           &in.composition[i][stress_start_index+n_independent_components]+n_independent_components));

                // $\eta^{t}_{effcreep}$. This viscosity has been calculated with the timestep_ratio dtc/dte.
                const double effective_creep_viscosity = effective_creep_viscosities[i];

                // $\eta_{el} = G \Delta t_c$.
                // The elastic viscosity has also already been scaled with $\frac{\Delta t_c} / {\Delta t_{el}}$
                // in light of the linear interpolation between $t$ and $t+ \Delta t_{el}$
                // when  $\Delta t_c$ and $t+\Delta t_el$ differ.
                const double elastic_viscosity = calculate_elastic_viscosity(average_elastic_shear_moduli[i]);

                // The ratio between the computational and elastic timestep $\frac{\Delta t_c} / {\Delta t_{el}}$.
                const double timestep_ratio = calculate_timestep_ratio();

                // Compute the total stress at time t.
                const SymmetricTensor<2, dim>
                stress_t = 2. * effective_creep_viscosity * Utilities::Tensors::consistent_deviator(in.strain_rate[i])
                           + effective_creep_viscosity / elastic_viscosity * stress_0_t
                           + (1. - timestep_ratio) * (1. - effective_creep_viscosity / elastic_viscosity) * stress_old;

                // Fill reaction rates.
                // During this timestep, the reaction rates will be multiplied
                // with the current timestep size to turn the rate of change into a change.
                // However, this update belongs
                // to the previous timestep. Therefore we divide by the
                // current timestep and multiply with the previous one.
                // When multiplied with the current timestep, this will give
                // (rate * previous_dt / current_dt) * current_dt = rate * previous_dt = previous_change.
                // previous_change = stress_t - stress_0_t.
                // To compute the rate we should return to the operator splitting scheme,
                // we therefore divide the change in stress by the current timestep current_dt (= dtc).
                const double dtc = timestep_ratio * elastic_timestep();

                const SymmetricTensor<2, dim> stress_update = (stress_t - stress_0_t) / dtc;

                Utilities::Tensors::unroll_symmetric_tensor_into_array(stress_update,
                                                                       &reaction_rate_out->reaction_rates[i][stress_start_index],
                                                                       &reaction_rate_out->reaction_rates[i][stress_start_index]+n_independent_components);

                // Also update the second set of stresses, stress_old, with the newly computed stress,
                // which in the rest of the timestep will serve as the old stress advected but not rotated
                // into the current timestep. This function fill_reaction_rates is only called at the
                // beginning of the timestep, and so this update only happens once.
                const SymmetricTensor<2, dim> stress_old_update = (stress_t - stress_old) / dtc;

                Utilities::Tensors::unroll_symmetric_tensor_into_array(stress_old_update,
                                                                       &reaction_rate_out->reaction_rates[i][stress_start_index+n_independent_components],
                                                                       &reaction_rate_out->reaction_rates[i][stress_start_index+n_independent_components]+n_independent_components);
              }
          }
      }



      template <int dim>
      double
      Elasticity<dim>::elastic_timestep () const
      {
        // The elastic time step ($\Delta t_el$, dte) is equal to the numerical time step if the time step number
        // is greater than 0 and the parameter 'use_fixed_elastic_time_step' is set to false.
        // On the first (0) time step, the elastic time step is always equal to the value
        // specified in 'fixed_elastic_time_step', which is also used in all subsequent time
        // steps if 'use_fixed_elastic_time_step' is set to true.
        //
        // We also use this parameter when we are still *before* the first time step,
        // i.e., if the time step number is numbers::invalid_unsigned_int.
        if (use_fixed_elastic_time_step && this->get_timestep_number() > 0 && this->simulator_is_past_initialization())
          AssertThrow(fixed_elastic_time_step >= this->get_timestep(), ExcMessage("The elastic timestep has to be equal to or bigger than the numerical timestep"));

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
      double
      Elasticity<dim>::calculate_timestep_ratio() const
      {
        // Before the simulator is initialized, get_timestep() can return
        // a NaN. During assembly in timestep 0, get_timestep() returns 0. Therefore we guess the
        // timestep size using the maximum timestep parameters capped by the elastic timestep.
        double dtc = this->get_timestep();
        if (!this->simulator_is_past_initialization() ||
            (this->get_timestep_number() == 0 && this->get_timestep() == 0))
          dtc = std::min(std::min(this->get_parameters().maximum_time_step, this->get_parameters().maximum_first_time_step), elastic_timestep());

        return dtc / elastic_timestep();
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
        // In case the computational timestep dtc differs from the elastic timestep,
        // we need to scale the elastic viscosity with the timestep ratio dtc/dte:
        // $\eta_{el} = \Delta t_{el} G$
        // $\eta_{el}^{c} = \Delta t_{el} G \frac{\Delta t_c}{\Delta t_{el}} = \Delta t_c G$.
        // Since we already have a function that returns the timestep ratio $\frac{\Delta t_c}{\Delta t_{el}}$
        // and deals with the special cases where $\Delta t_c$ is not know yet, use that function.
        // The damper is asserted to be zero when using this rheology at the moment. If the
        // damper were to be applied, the user should input the damper viscosity over
        // the computational timestep.
        const double timestep_ratio = calculate_timestep_ratio();
        return shear_modulus*elastic_timestep()*timestep_ratio + elastic_damper_viscosity;
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_viscoelastic_viscosity (const double viscosity,
                                        const double shear_modulus) const
      {
        // The elastic viscosity has been scaled with the timestep ratio $\frac{\Delta t_c}{\Delta t_{el}}$.
        // The viscous viscosity has not been scaled yet, so we do it here.
        // Scaling both viscosities with the timestep ratio before computing the effective
        // viscoelastic viscosity equals scaling the viscoelastic viscosity computed from
        // the unscaled elastic and viscous viscosity.
        // If the computational timestep is small compared to the elastic timestep,
        // the effective viscosity becomes small as well. This reduction is balanced by
        // an increasing body force term in the right-hand-side of the momentum equation.
        const double timestep_ratio = calculate_timestep_ratio();
        const double elastic_viscosity = calculate_elastic_viscosity(shear_modulus);
        return 1. / (1./elastic_viscosity + 1./(viscosity*timestep_ratio));
      }



      template <int dim>
      SymmetricTensor<2,dim>
      Elasticity<dim>::
      calculate_viscoelastic_strain_rate(const SymmetricTensor<2,dim> &strain_rate,
                                         const SymmetricTensor<2,dim> &stress_0_advected,
                                         const SymmetricTensor<2,dim> &stress_old,
                                         const double viscosity_pre_yield,
                                         const double shear_modulus) const
      {
        // The first term in the following expression is the overall strain rate
        // of one or more isostress rheological elements (in series).
        // One of these elements must be an elastic component (potentially damped).
        // The second term corresponds to a fictional strain rate arising from
        // elastic stresses stored from the last time step. Note the parallels with the
        // viscous part of the strain rate deviator,
        // which is equal to 0.5 * stress / viscosity.

        // The elastic viscosity is already scaled with the timestep ratio.
        const double elastic_viscosity = calculate_elastic_viscosity(shear_modulus);
        // viscosity_pre_yield is also already scaled with the timestep ratio.
        const double creep_viscosity = viscosity_pre_yield;

        // The ratio between the computational and elastic timestep.
        const double timestep_ratio = calculate_timestep_ratio();

        const SymmetricTensor<2, dim>
        edot = strain_rate + 0.5 * stress_0_advected / elastic_viscosity
               + 0.5 * (1. - timestep_ratio) * (1.  - creep_viscosity/elastic_viscosity) * stress_old / creep_viscosity;

        return edot;
      }



      template <int dim>
      std::vector<SymmetricTensor<2, dim>>
      Elasticity<dim>::
      retrieve_stress_previous_timestep (const MaterialModel::MaterialModelInputs<dim> &in,
                                         const std::vector<Point<dim>> &quadrature_positions) const
      {
        std::vector<SymmetricTensor<2, dim>> stress_t(in.n_evaluation_points(), SymmetricTensor<2, dim>());

        const unsigned int stress_start_index = this->introspection().compositional_index_for_name("ve_stress_xx");

        // Either get stress_t - the fully updated stress of the previous timestep - from the fields
        // or from the particles.
        if (this->get_parameters().compositional_field_methods[stress_start_index] == Parameters<dim>::AdvectionFieldMethod::fem_field)
          {
            // Get the compositional fields from the previous timestep $t$.
            // The 'old_solution' has been updated to the full stress tensor
            // of time $t$ by the operator splitting step at the beginning
            // of the current timestep.
            std::vector<double> old_solution_values(this->get_fe().dofs_per_cell);
            in.current_cell->get_dof_values(this->get_old_solution(),
                                            old_solution_values.begin(),
                                            old_solution_values.end());

            // Only create the evaluator the first time we get here.
            // With the field types, we can avoid specifying the stress
            // as the first fields. However, the FEPointEvaluation evaluators
            // can select a range of fields, but this has to be a consecutive
            // range starting from a given index. To avoid having to evaluate
            // all fields, we still request that the stress fields are listed
            // in a consecutive order without interruption by other fields.
            if (!evaluator_composition)
              evaluator_composition.reset(new FEPointEvaluation<n_independent_components, dim>(this->get_mapping(),
                                          this->get_fe(),
                                          update_values,
                                          this->introspection().component_indices.compositional_fields[stress_start_index]));

            // Initialize the evaluator for the composition values.
            evaluator_composition->reinit(in.current_cell, quadrature_positions);
            evaluator_composition->evaluate(old_solution_values,
                                            EvaluationFlags::values);

            // Get the composition values representing the viscoelastic stress field tensor components
            // of the previous timestep from the evaluator.
            // We assume (and assert in parse_parameters) that the 2*n_independent_components
            // tensor components are in the correct order and consecutively listed.
            // These stresses have not yet been rotated or advected to the current timestep.
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                const Tensor<1,n_independent_components> composition_values = evaluator_composition->get_value(i);
                stress_t[i] = Utilities::Tensors::to_symmetric_tensor<dim>(&composition_values[0],
                                                                           &composition_values[0]+n_independent_components);
              }
          }
        else
          {
            // If we use particles instead of fields to track the stresses, use the MaterialInputs,
            // which include the stresses stored on the particles.
            // The computation of the reaction terms (during which this function is called) happens
            // after the particles have been restored to their position and values
            // at the beginning of the timestep (i.e., the position and values of the previous timestep) and after they
            // have been updated to the full stress of the previous timestep by the reaction rates.
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              stress_t[i] = Utilities::Tensors::to_symmetric_tensor<dim>(&in.composition[i][stress_start_index],
                                                                         &in.composition[i][stress_start_index]+n_independent_components);
          }

        return stress_t;
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
