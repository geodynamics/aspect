/*
  Copyright (C) 2019 - 2020 by the authors of the ASPECT code.

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

#include <deal.II/fe/fe_values.h>
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
                           Patterns::List(Patterns::Double (0.)),
                           "List of elastic shear moduli, $G$, "
                           "for background material and compositional fields, "
                           "for a total of N+1 values, where N is the number of compositional fields. "
                           "The default value of 75 GPa is representative of mantle rocks. Units: Pa.");
        prm.declare_entry ("Use fixed elastic time step", "unspecified",
                           Patterns::Selection("true|false|unspecified"),
                           "Select whether the material time scale in the viscoelastic constitutive "
                           "relationship uses the regular numerical time step or a separate fixed "
                           "elastic time step throughout the model run. The fixed elastic time step "
                           "is always used during the initial time step. If a fixed elastic time "
                           "step is used throughout the model run, a stress averaging scheme can be "
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
        prm.declare_entry ("Use stress averaging","false",
                           Patterns::Bool (),
                           "Whether to apply a stress averaging scheme to account for differences "
                           "between the fixed elastic time step and numerical time step. ");
      }



      template <int dim>
      void
      Elasticity<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Get the number of fields for composition-dependent material properties
        const unsigned int n_fields = this->n_compositional_fields() + 1;

        elastic_shear_moduli = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Elastic shear moduli"))),
                                                                       n_fields,
                                                                       "Elastic shear moduli");

        if (prm.get ("Use fixed elastic time step") == "true")
          use_fixed_elastic_time_step = true;
        else if (prm.get ("Use fixed elastic time step") == "false")
          use_fixed_elastic_time_step = false;
        else
          AssertThrow(false, ExcMessage("'Use fixed elastic time step' must be set to 'true' or 'false'"));

        use_stress_averaging = prm.get_bool ("Use stress averaging");
        if (use_stress_averaging)
          AssertThrow(use_fixed_elastic_time_step == true,
                      ExcMessage("Stress averaging can only be used if 'Use fixed elastic time step' is set to true'"));

        fixed_elastic_time_step = prm.get_double ("Fixed elastic time step");
        AssertThrow(fixed_elastic_time_step > 0,
                    ExcMessage("The fixed elastic time step must be greater than zero"));

        if (this->convert_output_to_years())
          fixed_elastic_time_step *= year_in_seconds;

        AssertThrow(this->get_parameters().enable_elasticity == true,
                    ExcMessage ("Material model Viscoelastic only works if 'Enable elasticity' is set to true"));

        // Check whether the compositional fields representing the viscoelastic
        // stress tensor are both named correctly and listed in the right order.
        AssertThrow(this->introspection().compositional_index_for_name("stress_xx") == 0,
                    ExcMessage("Rheology model Elasticity only works if the first "
                               "compositional field is called stress_xx."));
        AssertThrow(this->introspection().compositional_index_for_name("stress_yy") == 1,
                    ExcMessage("Rheology model Elasticity only works if the second "
                               "compositional field is called stress_yy."));
        if (dim == 2)
          {
            AssertThrow(this->introspection().compositional_index_for_name("stress_xy") == 2,
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field is called stress_xy."));
          }
        else if (dim == 3)
          {
            AssertThrow(this->introspection().compositional_index_for_name("stress_zz") == 2,
                        ExcMessage("Rheology model Elasticity only works if the third "
                                   "compositional field is called stress_zz."));
            AssertThrow(this->introspection().compositional_index_for_name("stress_xy") == 3,
                        ExcMessage("Rheology model Elasticity only works if the fourth "
                                   "compositional field is called stress_xy."));
            AssertThrow(this->introspection().compositional_index_for_name("stress_xz") == 4,
                        ExcMessage("Rheology model Elasticity only works if the fifth "
                                   "compositional field is called stress_xz."));
            AssertThrow(this->introspection().compositional_index_for_name("stress_yz") == 5,
                        ExcMessage("Rheology model Elasticity only works if the sixth "
                                   "compositional field is called stress_yz."));
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
                     Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes),
                    ExcMessage("The material model will only work with the nonlinear "
                               "solver schemes 'single Advection, single Stokes' and "
                               "'single Advection, iterated Stokes'"));

        // Functionality to average the additional RHS terms over the cell is not implemented.
        // Consequently, it is only possible to use elasticity with the Material averaging schemes
        // 'none', 'harmonic average only viscosity', 'project to Q1 only viscosity'.
        AssertThrow((this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::none
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::harmonic_average_only_viscosity
                     ||
                     this->get_parameters().material_averaging == MaterialModel::MaterialAveraging::project_to_Q1_only_viscosity),
                    ExcMessage("Material models with elasticity can only be used with the material "
                               "averaging schemes 'none', 'harmonic average only viscosity', and "
                               "project to Q1 only viscosity'. This parameter ('Material averaging') "
                               "is located within the 'Material model' subsection."));
      }



      template <int dim>
      void
      Elasticity<dim>::create_elastic_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<ElasticAdditionalOutputs<dim> >() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<ElasticAdditionalOutputs<dim>> (n_points));
          }
      }



      template <int dim>
      void
      Elasticity<dim>::fill_elastic_force_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                   const std::vector<double> &average_elastic_shear_moduli,
                                                   MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // Create a reference to the structure for the elastic force terms that are needed to compute the
        // right-hand side of the Stokes system
        MaterialModel::ElasticOutputs<dim>
        *force_out = out.template get_additional_output<MaterialModel::ElasticOutputs<dim> >();

        if (force_out == nullptr)
          return;

        if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0 && in.requests_property(MaterialProperties::reaction_terms))
          {
            const double dte = elastic_timestep();

            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                // Get old stresses from compositional fields
                SymmetricTensor<2,dim> stress_old;
                for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                  stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];

                // Average viscoelastic viscosity
                const double average_viscoelastic_viscosity = out.viscosities[i];

                // Fill elastic force outputs (See equation 30 in Moresi et al., 2003, J. Comp. Phys.)
                force_out->elastic_force[i] = -1. * ( ( average_viscoelastic_viscosity / ( average_elastic_shear_moduli[i] * dte  ) ) * stress_old );

              }
          }
        else
          std::fill(force_out->elastic_force.begin(), force_out->elastic_force.end(), 0.0);
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
            std::vector<Point<dim> > quadrature_positions(in.n_evaluation_points());
            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

            // FEValues requires a quadrature and we provide the default quadrature
            // as we only need to evaluate the solution and gradients.
            FEValues<dim> fe_values (this->get_mapping(),
                                     this->get_fe(),
                                     Quadrature<dim>(quadrature_positions),
                                     update_gradients);

            fe_values.reinit (in.current_cell);
            std::vector<Tensor<2,dim> > old_velocity_gradients (quadrature_positions.size(), Tensor<2,dim>());
            fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_old_solution(),
                                                                                           old_velocity_gradients);

            const double dte = elastic_timestep();
            const double dt = this->get_timestep();

            for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
              {
                // Get old stresses from compositional fields
                SymmetricTensor<2,dim> stress_old;
                for (unsigned int j=0; j < SymmetricTensor<2,dim>::n_independent_components; ++j)
                  stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)] = in.composition[i][j];

                // Calculate the rotated stresses
                // Rotation (vorticity) tensor (equation 25 in Moresi et al., 2003, J. Comp. Phys.)
                const Tensor<2,dim> rotation = 0.5 * ( old_velocity_gradients[i] - transpose(old_velocity_gradients[i]) );

                // Average viscoelastic viscosity
                const double average_viscoelastic_viscosity = out.viscosities[i];

                // Calculate the current (new) viscoelastic stress, which is a function of the material
                // properties (viscoelastic viscosity, shear modulus), elastic time step size, strain rate,
                // vorticity and prior (inherited) viscoelastic stresses (see equation 29 in Moresi et al.,
                // 2003, J. Comp. Phys.)
                SymmetricTensor<2,dim> stress_new = ( 2. * average_viscoelastic_viscosity * deviator(in.strain_rate[i]) ) +
                                                    ( ( average_viscoelastic_viscosity / ( average_elastic_shear_moduli[i] * dte ) ) * stress_old ) +
                                                    ( ( average_viscoelastic_viscosity / average_elastic_shear_moduli[i] ) *
                                                      ( symmetrize(rotation * Tensor<2,dim>(stress_old) ) - symmetrize(Tensor<2,dim>(stress_old) * rotation) ) );

                // Stress averaging scheme to account for difference between fixed elastic time step
                // and numerical time step (see equation 32 in Moresi et al., 2003, J. Comp. Phys.)
                if (use_fixed_elastic_time_step == true && use_stress_averaging == true)
                  {
                    stress_new = ( ( 1. - ( dt / dte ) ) * stress_old ) + ( ( dt / dte ) * stress_new ) ;
                  }

                // Fill reaction terms
                for (unsigned int j = 0; j < SymmetricTensor<2,dim>::n_independent_components ; ++j)
                  out.reaction_terms[i][j] = -stress_old[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)]
                                             + stress_new[SymmetricTensor<2,dim>::unrolled_to_component_indices(j)];

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
                             this->get_timestep()
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
      calculate_viscoelastic_viscosity (const double viscosity,
                                        const double elastic_shear_modulus) const
      {
        const double elastic_viscosity = elastic_shear_modulus*elastic_timestep();
        return 1. / (1./elastic_viscosity + 1./viscosity);
      }



      template <int dim>
      double
      Elasticity<dim>::
      calculate_viscoelastic_strain_rate(const SymmetricTensor<2,dim> &strain_rate,
                                         const SymmetricTensor<2,dim> &stress,
                                         const double shear_modulus) const
      {
        // The second term in the following expression corresponds to the
        // elastic part of the strain rate deviator. Note the parallels with the
        // viscous part of the strain rate deviator,
        // which is equal to 0.5 * stress / viscosity.
        const SymmetricTensor<2,dim> edot_deviator = deviator(strain_rate) + 0.5*stress /
                                                     (shear_modulus * elastic_timestep());

        return std::sqrt(std::fabs(second_invariant(edot_deviator)));
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
