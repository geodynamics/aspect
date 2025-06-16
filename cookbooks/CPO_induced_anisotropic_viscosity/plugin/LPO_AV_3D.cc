
/*
 Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include "LPO_AV_3D.h"
#include <aspect/material_model/simple.h>
#include <aspect/material_model/grain_size.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/heating_model/shear_heating.h>
#include <aspect/heating_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator_signals.h>
#include <aspect/postprocess/particles.h>
#include <aspect/introspection.h>
#include <aspect/plugins.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/lapack_templates.h>
#include <deal.II/lac/scalapack.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/physics/notation.h>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS


namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for anisotropic viscosities to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */



    namespace
    {

      template <int dim>
      std::vector<std::string> make_AV_additional_outputs_names()
      {
        std::vector<std::string> names;

        for (unsigned int i = 0; i < Tensor<4,dim>::n_independent_components ; ++i)
          {
            TableIndices<4> indices(Tensor<4,dim>::unrolled_to_component_indices(i));
            names.push_back("anisotropic_viscosity"+std::to_string(indices[0]+1)+std::to_string(indices[1]+1)+std::to_string(indices[2]+1)+std::to_string(indices[3]+1));
          }
        return names;
      }
    }



    template <int dim>
    AV<dim>::AV (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_AV_additional_outputs_names<dim>()),
      stress_strain_directors(n_points, dealii::identity_tensor<dim> ())
    {}



    template <int dim>
    std::vector<double>
    AV<dim>::get_nth_output(const unsigned int idx) const
    {
      std::vector<double> output(stress_strain_directors.size());
      for (unsigned int i = 0; i < stress_strain_directors.size() ; ++i)
        {
          output[i]= stress_strain_directors[i][Tensor<4,dim>::unrolled_to_component_indices(idx)];
        }
      return output;
    }
  }
}

namespace aspect
{
  namespace Assemblers
  {
    /**
     * A class containing the functions to assemble the Stokes preconditioner.
     */
    template <int dim>
    class StokesPreconditionerAV : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const override;

        /**
         * Create AnisotropicViscosities.
         */
        void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Stokes equation for the current cell.
     */
    template <int dim>
    class StokesIncompressibleTermsAV : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const override;

        /**
         * Create AdditionalMaterialOutputsStokesRHS if we need to do so.
         */
        void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const override;
    };



    template <int dim>
    void
    StokesPreconditionerAV<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>&> (data_base);

      const std::shared_ptr<const MaterialModel::AV<dim>> anisotropic_viscosity
        = scratch.material_model_outputs.template get_additional_output_object<MaterialModel::AV<dim>>();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points           = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      // First loop over all dofs and find those that are in the Stokes system
      // save the component (pressure and dim velocities) each belongs to.
      for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
        {
          if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
            {
              scratch.dof_component_indices[i_stokes] = fe.system_to_component_index(i).first;
              ++i_stokes;
            }
          ++i;
        }

      // Loop over all quadrature points and assemble their contributions to
      // the preconditioner matrix
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          for (unsigned int i = 0, i_stokes = 0; i_stokes < stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.grads_phi_u[i_stokes] =
                    scratch.finite_element_values[introspection.extractors
                                                  .velocities].symmetric_gradient(i, q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection
                                                                          .extractors.pressure].value(i, q);
                  ++i_stokes;
                }
              ++i;
            }

          const double eta = scratch.material_model_outputs.viscosities[q];
          const double one_over_eta = 1. / eta;
          const SymmetricTensor<4, dim> &stress_strain_director = anisotropic_viscosity->stress_strain_directors[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += (2.0 * eta * (scratch.grads_phi_u[i]
                                                         * stress_strain_director
                                                         * scratch.grads_phi_u[j])
                                            + one_over_eta * pressure_scaling
                                            * pressure_scaling
                                            * (scratch.phi_p[i]
                                               * scratch.phi_p[j]))
                                           * JxW;
        }
    }



    template <int dim>
    void
    StokesPreconditionerAV<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template has_additional_output_object<MaterialModel::AV<dim>>() == false)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::AV<dim>> (n_points));
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAV<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>&> (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>&> (data_base);

      const std::shared_ptr<const MaterialModel::AV<dim>> anisotropic_viscosity
        = scratch.material_model_outputs.template get_additional_output_object<MaterialModel::AV<dim>>();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      const std::shared_ptr<const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> force
        = scratch.material_model_outputs.template get_additional_output_object<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>();

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].value (i,q);
                  scratch.phi_p[i_stokes] = scratch.finite_element_values[introspection.extractors.pressure].value (i, q);
                  if (scratch.rebuild_stokes_matrix)
                    {
                      scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);
                    }
                  ++i_stokes;
                }
              ++i;
            }
          // Viscosity scalar
          const double eta = (scratch.rebuild_stokes_matrix
                              ?
                              scratch.material_model_outputs.viscosities[q]
                              :
                              numbers::signaling_nan<double>());

          const SymmetricTensor<4, dim> &stress_strain_director = anisotropic_viscosity->stress_strain_directors[q];

          const Tensor<1,dim>
          gravity = this->get_gravity_model().gravity_vector (scratch.finite_element_values.quadrature_point(q));

          const double density = scratch.material_model_outputs.densities[q];
          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            {
              data.local_rhs(i) += (density * gravity * scratch.phi_u[i])
                                   * JxW;

              if (force != nullptr)
                data.local_rhs(i) += (force->rhs_u[q] * scratch.phi_u[i]
                                      + pressure_scaling * force->rhs_p[q] * scratch.phi_p[i])
                                     * JxW;

              if (scratch.rebuild_stokes_matrix)
                for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
                  {
                    data.local_matrix(i,j) += ( eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                                // assemble \nabla p as -(p, div v):
                                                - (pressure_scaling *
                                                   scratch.div_phi_u[i] * scratch.phi_p[j])
                                                // assemble the term -div(u) as -(div u, q).
                                                // Note the negative sign to make this
                                                // operator adjoint to the grad p term:
                                                - (pressure_scaling *
                                                   scratch.phi_p[i] * scratch.div_phi_u[j]))
                                              * JxW;
                  }
            }
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAV<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template has_additional_output_object<MaterialModel::AV<dim>>() == false)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::AV<dim>> (n_points));
        }

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template has_additional_output_object<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>() == false)
        {
          outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));
        }
      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output_object<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>>()->rhs_u.size()
             == n_points, ExcInternalError());
    }
  }

  namespace HeatingModel
  {
    template <int dim>
    class ShearHeatingAV : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        /**
         * Allow the heating model to attach additional material model outputs.
         */
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const override;
    };



    template <int dim>
    void
    ShearHeatingAV<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      // Some material models provide dislocation viscosities and boundary area work fractions
      // as additional material outputs. If they are attached, use them.
      const std::shared_ptr<const ShearHeatingOutputs<dim>> shear_heating_out
        = material_model_outputs.template get_additional_output_object<ShearHeatingOutputs<dim>>();

      const std::shared_ptr<const MaterialModel::AV<dim>> anisotropic_viscosity
        = material_model_outputs.template get_additional_output_object<MaterialModel::AV<dim>>();

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // If there is an anisotropic viscosity, use it to compute the correct stress
          const SymmetricTensor<2,dim> &directed_strain_rate = ((anisotropic_viscosity != nullptr)
                                                                ?
                                                                anisotropic_viscosity->stress_strain_directors[q]
                                                                * material_model_inputs.strain_rate[q]
                                                                :
                                                                material_model_inputs.strain_rate[q]);

          const SymmetricTensor<2,dim> stress =
            2 * material_model_outputs.viscosities[q] *
            (this->get_material_model().is_compressible()
             ?
             directed_strain_rate - 1./3. * trace(directed_strain_rate) * unit_symmetric_tensor<dim>()
             :
             directed_strain_rate);

          const SymmetricTensor<2,dim> deviatoric_strain_rate =
            (this->get_material_model().is_compressible()
             ?
             material_model_inputs.strain_rate[q]
             - 1./3. * trace(material_model_inputs.strain_rate[q]) * unit_symmetric_tensor<dim>()
             :
             material_model_inputs.strain_rate[q]);

          heating_model_outputs.heating_source_terms[q] = stress * deviatoric_strain_rate;

          // If shear heating work fractions are provided, reduce the
          // overall heating by this amount (which is assumed to be converted into other forms of energy)
          if (shear_heating_out != nullptr)
            heating_model_outputs.heating_source_terms[q] *= shear_heating_out->shear_heating_work_fractions[q];

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ShearHeatingAV<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      const unsigned int n_points = material_model_outputs.viscosities.size();

      if (material_model_outputs.template has_additional_output_object<MaterialModel::AV<dim>>() == false)
        {
          material_model_outputs.additional_outputs.push_back(
            std::make_unique<MaterialModel::AV<dim>> (n_points));
        }

      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }
  }

}

namespace aspect
{

//Next session is a more evolved implementation of anisotropic viscosity in the material model based on Hansen et al 2016 and Kiraly et al 2020
  namespace MaterialModel
  {

    template <int dim>
    void
    LPO_AV_3D<dim>::set_assemblers(const SimulatorAccess<dim> &,
                                   Assemblers::Manager<dim> &assemblers) const
    {
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesPreconditioner<dim>>(*(assemblers.stokes_preconditioner[i])))
            assemblers.stokes_preconditioner[i] = std::make_unique<Assemblers::StokesPreconditionerAV<dim>> ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesIncompressibleTerms<dim>>(*(assemblers.stokes_system[i])))
            assemblers.stokes_system[i] = std::make_unique<Assemblers::StokesIncompressibleTermsAV<dim>> ();
        }
    }



    template <int dim>
    void
    LPO_AV_3D<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&LPO_AV_3D<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
      AssertThrow((dim==3),
                  ExcMessage("Olivine has 3 independent slip systems, allowing for deformation in 3 independent directions, hence these models only work in 3D"));

      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("phi1"));
      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("eigvalue_a1"));
      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("eigvalue_a2"));
      cpo_bingham_avg_a.push_back (this->introspection().compositional_index_for_name("eigvalue_a3"));

      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("theta"));
      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("eigvalue_b1"));
      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("eigvalue_b2"));
      cpo_bingham_avg_b.push_back (this->introspection().compositional_index_for_name("eigvalue_b3"));

      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("phi2"));
      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("eigvalue_c1"));
      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("eigvalue_c2"));
      cpo_bingham_avg_c.push_back (this->introspection().compositional_index_for_name("eigvalue_c3"));


    }



    template<int dim>
    Tensor<2,3>
    AV<dim>::euler_angles_to_rotation_matrix(double phi1, double theta, double phi2)
    {
      Tensor<2,3> rot_matrix;
      //R3*R2*R1 ZXZ rotation. Note it is not exactly the same as in utilities.cc
      rot_matrix[0][0] = cos(phi2)*cos(phi1) - cos(theta)*sin(phi1)*sin(phi2); //
      rot_matrix[0][1] = -cos(phi2)*sin(phi1) - cos(theta)*cos(phi1)*sin(phi2); //cos(phi2)*sin(phi1) + cos(theta)*cos(phi1)*sin(phi2);
      rot_matrix[0][2] = sin(phi2)*sin(theta);
      rot_matrix[1][0] = sin(phi2)*cos(phi1) + cos(theta)*sin(phi1)*cos(phi2); //-sin(phi2)*cos(phi1) - cos(theta)*sin(phi1)*cos(phi2);
      rot_matrix[1][1] = -sin(phi2)*sin(phi1) + cos(theta)*cos(phi1)*cos(phi2);
      rot_matrix[1][2] = -cos(phi2)*sin(theta); //cos(phi2)*sin(theta);
      rot_matrix[2][0] = sin(theta)*sin(phi1);
      rot_matrix[2][1] = sin(theta)*cos(phi1); //-sin(theta)*cos(phi1);
      rot_matrix[2][2] = cos(theta); //
      AssertThrow(rot_matrix[2][2] <= 1.0, ExcMessage("rot_matrix[2][2] > 1.0"));
      return rot_matrix;
    }



    template <>
    void
    LPO_AV_3D<2>::evaluate (const MaterialModel::MaterialModelInputs<2> &,
                            MaterialModel::MaterialModelOutputs<2> &) const
    {
      Assert (false, ExcNotImplemented());
    }


    template <>
    void
    LPO_AV_3D<3>::evaluate (const MaterialModel::MaterialModelInputs<3> &in,
                            MaterialModel::MaterialModelOutputs<3> &out) const
    {
      const int dim=3;
      const double sqrt2 = sqrt(2);
      MaterialModel::AV<dim> *anisotropic_viscosity
        = out.template get_additional_output<MaterialModel::AV<dim>>();
      EquationOfStateOutputs<dim> eos_outputs (1);



      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          // std::cout << "Evaluation point: " << q << std::endl;
          //change these according to diffusion dislocation material model I guess
          equation_of_state.evaluate(in, q, eos_outputs);

          // Get parameters for compute the effective viscosity
          // const double temperature = in.temperature[q];
          // const double pressure= in.pressure[q];
          const std::vector<double> composition = in.composition[q];
          const std::vector<double> volume_fractions = MaterialUtilities::compute_only_composition_fractions(composition,
                                                       this->introspection().chemical_composition_field_indices());

          // Compositional dependence of density, for Stokes sinker test
          const unsigned int comp_density_ind = this->introspection().compositional_index_for_name("D");
          double composition_dependence = composition[comp_density_ind];
          out.densities[q] = eos_outputs.densities[0] + compositional_delta_rho * composition_dependence;
          out.viscosities[q] = eta;
          out.thermal_expansion_coefficients[q] = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[q] = eos_outputs.specific_heat_capacities[0];
          out.thermal_conductivities[q] = 1;
          out.compressibilities[q] = eos_outputs.compressibilities[0];
          out.entropy_derivative_pressure[q] = eos_outputs.entropy_derivative_pressure[0];
          out.entropy_derivative_temperature[q] = eos_outputs.entropy_derivative_temperature[0];

          //Calculate effective viscosity
          const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
          const SymmetricTensor<2,dim> deviatoric_strain_rate
            = (this->get_material_model().is_compressible()
               ?
               strain_rate - 1./3. * trace(strain_rate) * unit_symmetric_tensor<dim>()
               :
               strain_rate);
          // std::cout << "deviatoric_strain_rate: " << deviatoric_strain_rate << std::endl;
          // The computation of the viscosity tensor is only necessary after the simulator has been initialized
          // and when the condition allows dislocation creep
          if  ((this->simulator_is_past_initialization()) && (this->get_timestep_number() > 0))
            {
              // std::cout << "T: " << in.temperature[q] << " det dev_sr: " << determinant(deviatoric_strain_rate) << std::endl;
              if ((in.temperature[q]>1000) && (isfinite(determinant(deviatoric_strain_rate))))// && (determinant(deviatoric_strain_rate) != 0)) && (isfinite(determinant(deviatoric_strain_rate)))
                {
                  // std::cout << "T: " << in.temperature[q] << " det dev_sr: " << determinant(deviatoric_strain_rate) << std::endl;
                  const unsigned int ind_vis = this->introspection().compositional_index_for_name("scalar_vis");
                  // std::cout << "Initial viscosity: " << composition[ind_vis] << std::endl;

                  //Create constant value to use for AV
                  const double A_o = 1.1e5*exp(-530000/(8.314*in.temperature[q]));
                  const double n = 3.5;
                  const double Gamma = (A_o/(std::pow(grain_size/1e6,0.73)));// in MPa^(-n)
                  // std::cout << "Gamma: " << Gamma << std::endl;
                  //SymmetricTensor<4,dim> old_stress_strain_director;
                  //std::vector<double> ssd_array(SymmetricTensor<4,dim>::n_independent_components);
                  //for (unsigned int i = 0; i < SymmetricTensor<4,dim>::n_independent_components ; ++i)
                  //{
                  //const unsigned int ind = this->introspection().compositional_index_for_name(ssd_names[i]);
                  //ssd_array[i] = composition[ind];
                  //AssertThrow(isfinite(composition[ind]),
                  //ExcMessage("Assigned prescribed field should be finite"));
                  //}
                  //std::copy(ssd_array.begin(), ssd_array.end(), old_stress_strain_director.begin_raw());

                  //Get eigen values from compositional fields
                  const double eigvalue_a1 = composition[cpo_bingham_avg_a[1]];
                  const double eigvalue_b1 = composition[cpo_bingham_avg_b[1]];
                  const double eigvalue_c1 = composition[cpo_bingham_avg_c[1]];
                  const double eigvalue_a2 = composition[cpo_bingham_avg_a[2]];
                  const double eigvalue_b2 = composition[cpo_bingham_avg_b[2]];
                  const double eigvalue_c2 = composition[cpo_bingham_avg_c[2]];
                  const double eigvalue_a3 = composition[cpo_bingham_avg_a[3]];
                  const double eigvalue_b3 = composition[cpo_bingham_avg_b[3]];
                  const double eigvalue_c3 = composition[cpo_bingham_avg_c[3]];

                  //Calculate the rotation matrix from the euler angles
                  const double phi1 = composition[cpo_bingham_avg_a[0]];
                  const double theta = composition[cpo_bingham_avg_b[0]];
                  const double phi2 = composition[cpo_bingham_avg_c[0]];
                  // std::cout<<"in mm: phi1 "<<phi1<<" theta "<<theta<<" phi2 "<<phi2<<std::endl;
                  Tensor<2,3> R = transpose(AV<dim>::euler_angles_to_rotation_matrix(phi1, theta, phi2));
                  // Tensor<2,3> R = transpose(Utilities::zxz_euler_angles_to_rotation_matrix(phi1*constants::radians_to_degree, theta*constants::radians_to_degree, phi2*constants::radians_to_degree));
                  // std::cout << "R in mm " << R << std::endl;

                  // // Check if the rotation matrix is orthogonal using R^T=R^(-1) and det(R)=1
                  // std::cout<<"in mm: transpose(R)*R= "<<transpose(R)*R<<std::endl;
                  // std::cout<<"in mm: transpose(R) "<<transpose(R)<<" invert(R) "<<invert(R)<<std::endl;
                  // std::cout<<"in mm: det(R)==1 "<<std::endl;

                  //Compute Hill Parameters FGHLMN from the eigenvalues of a,b,c axis
                  double F, G, H, L, M, N;
                  // F = std::pow(eigvalue_a1,2)*CnI_F[0] + eigvalue_a1*CnI_F[1] + eigvalue_a2*CnI_F[2] + (1/eigvalue_a3)*CnI_F[3] + std::pow(eigvalue_b1,2)*CnI_F[4] + eigvalue_b1*CnI_F[5] + eigvalue_b2*CnI_F[6] + (1/eigvalue_b3)*CnI_F[7] + std::pow(eigvalue_c1,2)*CnI_F[8] + eigvalue_c1*CnI_F[9] + eigvalue_c2*CnI_F[10] + (1/eigvalue_c3)*CnI_F[11] + CnI_F[12];
                  // G = std::pow(eigvalue_a1,2)*CnI_G[0] + eigvalue_a1*CnI_G[1] + eigvalue_a2*CnI_G[2] + (1/eigvalue_a3)*CnI_G[3] + std::pow(eigvalue_b1,2)*CnI_G[4] + eigvalue_b1*CnI_G[5] + eigvalue_b2*CnI_G[6] + (1/eigvalue_b3)*CnI_G[7] + std::pow(eigvalue_c1,2)*CnI_G[8] + eigvalue_c1*CnI_G[9] + eigvalue_c2*CnI_G[10] + (1/eigvalue_c3)*CnI_G[11] + CnI_G[12];
                  // H = std::pow(eigvalue_a1,2)*CnI_H[0] + eigvalue_a1*CnI_H[1] + eigvalue_a2*CnI_H[2] + (1/eigvalue_a3)*CnI_H[3] + std::pow(eigvalue_b1,2)*CnI_H[4] + eigvalue_b1*CnI_H[5] + eigvalue_b2*CnI_H[6] + (1/eigvalue_b3)*CnI_H[7] + std::pow(eigvalue_c1,2)*CnI_H[8] + eigvalue_c1*CnI_H[9] + eigvalue_c2*CnI_H[10] + (1/eigvalue_c3)*CnI_H[11] + CnI_H[12];
                  // L = std::abs(std::pow(eigvalue_a1,2)*CnI_L[0] + eigvalue_a1*CnI_L[1] + eigvalue_a2*CnI_L[2] + (1/eigvalue_a3)*CnI_L[3] + std::pow(eigvalue_b1,2)*CnI_L[4] + eigvalue_b1*CnI_L[5] + eigvalue_b2*CnI_L[6] + (1/eigvalue_b3)*CnI_L[7] + std::pow(eigvalue_c1,2)*CnI_L[8] + eigvalue_c1*CnI_L[9] + eigvalue_c2*CnI_L[10] + (1/eigvalue_c3)*CnI_L[11] + CnI_L[12]);
                  // M = std::abs(std::pow(eigvalue_a1,2)*CnI_M[0] + eigvalue_a1*CnI_M[1] + eigvalue_a2*CnI_M[2] + (1/eigvalue_a3)*CnI_M[3] + std::pow(eigvalue_b1,2)*CnI_M[4] + eigvalue_b1*CnI_M[5] + eigvalue_b2*CnI_M[6] + (1/eigvalue_b3)*CnI_M[7] + std::pow(eigvalue_c1,2)*CnI_M[8] + eigvalue_c1*CnI_M[9] + eigvalue_c2*CnI_M[10] + (1/eigvalue_c3)*CnI_M[11] + CnI_M[12]);
                  // N = std::abs(std::pow(eigvalue_a1,2)*CnI_N[0] + eigvalue_a1*CnI_N[1] + eigvalue_a2*CnI_N[2] + (1/eigvalue_a3)*CnI_N[3] + std::pow(eigvalue_b1,2)*CnI_N[4] + eigvalue_b1*CnI_N[5] + eigvalue_b2*CnI_N[6] + (1/eigvalue_b3)*CnI_N[7] + std::pow(eigvalue_c1,2)*CnI_N[8] + eigvalue_c1*CnI_N[9] + eigvalue_c2*CnI_N[10] + (1/eigvalue_c3)*CnI_N[11] + CnI_N[12]);

                  // old
                  F = std::pow(eigvalue_a1,2)*CnI_F[0] + eigvalue_a2*CnI_F[1] + (1/eigvalue_a3)*CnI_F[2] + std::pow(eigvalue_b1,2)*CnI_F[3] + eigvalue_b2*CnI_F[4] + (1/eigvalue_b3)*CnI_F[5] + std::pow(eigvalue_c1,2)*CnI_F[6] + eigvalue_c2*CnI_F[7] + (1/eigvalue_c3)*CnI_F[8] + CnI_F[9];
                  G = std::pow(eigvalue_a1,2)*CnI_G[0] + eigvalue_a2*CnI_G[1] + (1/eigvalue_a3)*CnI_G[2] + std::pow(eigvalue_b1,2)*CnI_G[3] + eigvalue_b2*CnI_G[4] + (1/eigvalue_b3)*CnI_G[5] + std::pow(eigvalue_c1,2)*CnI_G[6] + eigvalue_c2*CnI_G[7] + (1/eigvalue_c3)*CnI_G[8] + CnI_G[9];
                  H = std::pow(eigvalue_a1,2)*CnI_H[0] + eigvalue_a2*CnI_H[1] + (1/eigvalue_a3)*CnI_H[2] + std::pow(eigvalue_b1,2)*CnI_H[3] + eigvalue_b2*CnI_H[4] + (1/eigvalue_b3)*CnI_H[5] + std::pow(eigvalue_c1,2)*CnI_H[6] + eigvalue_c2*CnI_H[7] + (1/eigvalue_c3)*CnI_H[8] + CnI_H[9];
                  L = std::abs(std::pow(eigvalue_a1,2)*CnI_L[0] + eigvalue_a2*CnI_L[1] + (1/eigvalue_a3)*CnI_L[2] + std::pow(eigvalue_b1,2)*CnI_L[3] + eigvalue_b2*CnI_L[4] + (1/eigvalue_b3)*CnI_L[5] + std::pow(eigvalue_c1,2)*CnI_L[6] + eigvalue_c2*CnI_L[7] + (1/eigvalue_c3)*CnI_L[8] + CnI_L[9]);
                  M = std::abs(std::pow(eigvalue_a1,2)*CnI_M[0] + eigvalue_a2*CnI_M[1] + (1/eigvalue_a3)*CnI_M[2] + std::pow(eigvalue_b1,2)*CnI_M[3] + eigvalue_b2*CnI_M[4] + (1/eigvalue_b3)*CnI_M[5] + std::pow(eigvalue_c1,2)*CnI_M[6] + eigvalue_c2*CnI_M[7] + (1/eigvalue_c3)*CnI_M[8] + CnI_M[9]);
                  N = std::abs(std::pow(eigvalue_a1,2)*CnI_N[0] + eigvalue_a2*CnI_N[1] + (1/eigvalue_a3)*CnI_N[2] + std::pow(eigvalue_b1,2)*CnI_N[3] + eigvalue_b2*CnI_N[4] + (1/eigvalue_b3)*CnI_N[5] + std::pow(eigvalue_c1,2)*CnI_N[6] + eigvalue_c2*CnI_N[7] + (1/eigvalue_c3)*CnI_N[8] + CnI_N[9]);

                  // std::cout<<"in mm: eigvalue_a1 "<<eigvalue_a1<<" eigvalue_a2 "<<eigvalue_a2<<" eigvalue_a3 "<<eigvalue_a3<<std::endl;
                  // std::cout<<"mm: eigvalue_b1 "<<eigvalue_b1<<" eigvalue_b2 "<<eigvalue_b2<<" eigvalue_b3 "<<eigvalue_b3<<std::endl;
                  // std::cout<<"mm: eigvalue_c1 "<<eigvalue_c1<<" eigvalue_c2 "<<eigvalue_c2<<" eigvalue_c3 "<<eigvalue_c3<<std::endl;
                  // std::cout<<"F "<<F<<" G "<<G<<" H "<<H<<" L "<<L<<" M "<<M<<" N "<<N<<std::endl;
                  // F=0.5; G=0.5, H=0.5; L=1.5; M=1.5; N=1.5;
                  Tensor<2,6> R_CPO_K;
                  R_CPO_K[0][0] = std::pow(R[0][0],2);
                  R_CPO_K[0][1] = std::pow(R[0][1],2);
                  R_CPO_K[0][2] = std::pow(R[0][2],2);
                  R_CPO_K[0][3] = sqrt2*R[0][1]*R[0][2];
                  R_CPO_K[0][4] = sqrt2*R[0][0]*R[0][2];
                  R_CPO_K[0][5] = sqrt2*R[0][0]*R[0][1];

                  R_CPO_K[1][0] = std::pow(R[1][0],2);
                  R_CPO_K[1][1] = std::pow(R[1][1],2);
                  R_CPO_K[1][2] = std::pow(R[1][2],2);
                  R_CPO_K[1][3] = sqrt2*R[1][1]*R[1][2];
                  R_CPO_K[1][4] = sqrt2*R[1][0]*R[1][2];
                  R_CPO_K[1][5] = sqrt2*R[1][0]*R[1][1];

                  R_CPO_K[2][0] = std::pow(R[2][0],2);
                  R_CPO_K[2][1] = std::pow(R[2][1],2);
                  R_CPO_K[2][2] = std::pow(R[2][2],2);
                  R_CPO_K[2][3] = sqrt2*R[2][1]*R[2][2];
                  R_CPO_K[2][4] = sqrt2*R[2][0]*R[2][2];
                  R_CPO_K[2][5] = sqrt2*R[2][0]*R[2][1];

                  R_CPO_K[3][0] = sqrt2*R[1][0]*R[2][0];
                  R_CPO_K[3][1] = sqrt2*R[1][1]*R[2][1];
                  R_CPO_K[3][2] = sqrt2*R[1][2]*R[2][2];
                  R_CPO_K[3][3] = R[1][1]*R[2][2]+R[1][2]*R[2][1];
                  R_CPO_K[3][4] = R[1][0]*R[2][2]+R[1][2]*R[2][0];
                  R_CPO_K[3][5] = R[1][0]*R[2][1]+R[1][1]*R[2][0];

                  R_CPO_K[4][0] = sqrt2*R[0][0]*R[2][0];
                  R_CPO_K[4][1] = sqrt2*R[0][1]*R[2][1];
                  R_CPO_K[4][2] = sqrt2*R[0][2]*R[2][2];
                  R_CPO_K[4][3] = R[0][1]*R[2][2]+R[0][2]*R[2][1];
                  R_CPO_K[4][4] = R[0][0]*R[2][2]+R[0][2]*R[2][0];
                  R_CPO_K[4][5] = R[0][0]*R[2][1]+R[0][1]*R[2][0];

                  R_CPO_K[5][0] = sqrt2*R[0][0]*R[1][0];
                  R_CPO_K[5][1] = sqrt2*R[0][1]*R[1][1];
                  R_CPO_K[5][2] = sqrt2*R[0][2]*R[1][2];
                  R_CPO_K[5][3] = R[0][1]*R[1][2]+R[0][2]*R[1][1];
                  R_CPO_K[5][4] = R[0][0]*R[1][2]+R[0][2]*R[1][0];
                  R_CPO_K[5][5] = R[0][0]*R[1][1]+R[0][1]*R[1][0];

                  // // Check if the rotation matrix is orthogonal using R^T=R^(-1) and det(R)=1
                  // std::cout<<"in mm: transpose(R_CPO_K)*R_CPO_K= "<<transpose(R_CPO_K)*R_CPO_K<<std::endl;
                  // std::cout<<"in mm: transpose(R_CPO_K) "<<transpose(R_CPO_K)<<" invert(R_CPO_K) "<<invert(R)<<std::endl;
                  // std::cout<<"in mm: det(R)==1 "<<std::endl;

                  SymmetricTensor<2,6> A;
                  A[0][0] = 2.0/3.0*(F+H);
                  A[0][1] = 2.0/3.0*(-F);
                  A[0][2] = 2.0/3.0*(-H);
                  A[1][1] = 2.0/3.0*(G+F);
                  A[1][2] = 2.0/3.0*(-G);
                  A[2][2] = 2.0/3.0*(H+G);
                  A[3][3] = 2.0/3.0*L;
                  A[4][4] = 2.0/3.0*M;
                  A[5][5] = 2.0/3.0*N;
                  // std::cout << "A (fluidity) " << A << std::endl;

                  //Invert using ScaLAPACK in dealii
                  FullMatrix<double> A_mat(6, 6);
                  for (unsigned int ai=0; ai<6; ++ai)
                    {
                      for (unsigned int aj=0; aj<6; ++aj)
                        {
                          A_mat(ai,aj) = A[ai][aj];
                        }
                    }
                  const double ratio = 1e-8;
                  std::shared_ptr<Utilities::MPI::ProcessGrid> grid = std::make_shared<Utilities::MPI::ProcessGrid>(this->get_mpi_communicator(),6,6,4,4);
                  ScaLAPACKMatrix<double> A_scalapack(6,6,grid,4,4);
                  A_scalapack = A_mat;
                  A_scalapack.pseudoinverse(ratio);
                  FullMatrix<double> pinvA_mat(6,6);
                  A_scalapack.copy_to(pinvA_mat);

                  SymmetricTensor<2,6> invA;
                  for (unsigned int ai=0; ai<6; ++ai)
                    {
                      for (unsigned int aj=0; aj<6; ++aj)
                        {
                          invA[ai][aj] = pinvA_mat(ai,aj);
                        }
                    }
                  // std::cout << "invA " << invA <<std::endl;

                  //Calculate the fluidity tensor in the LPO frame
                  Tensor<2,6> V = R_CPO_K * invA * transpose(R_CPO_K);//invA;//
                  // Tensor<2,6> V = transpose(R_CPO_K) * invA * R_CPO_K;//invA;//
                  // std::cout << "V: " << V <<std::endl;

                  // SymmetricTensor<2,6> V;
                  // V[0][0] = 2.0/3.0;
                  // V[0][1] = -1.0/3.0;
                  // V[0][2] = -1.0/3.0;
                  // V[1][1] = 2.0/3.0;
                  // V[1][2] = -1.0/3.0;
                  // V[2][2] = 2.0/3.0;
                  // V[3][3] = 1;
                  // V[4][4] = 1;
                  // V[5][5] = 1;

                  //Convert rank 2 viscosity tensor to rank 4
                  FullMatrix<double> V_mat(6,6);
                  for (unsigned int vi=0; vi<6; ++vi)
                    {
                      for (unsigned int vj=0; vj<6; ++vj)
                        {
                          V_mat[vi][vj] = V[vi][vj];
                        }
                    }
                  SymmetricTensor<4,dim> V_r4;
                  dealii::Physics::Notation::Kelvin::to_tensor(V_mat, V_r4);

                  if (anisotropic_viscosity != nullptr)
                    {
                      anisotropic_viscosity->stress_strain_directors[q] = V_r4;
                    }

                  double scalar_viscosity = composition[ind_vis];

                  //In first time step using input viscosity can lead to convergence issue if the strainrate varies significantly within the model domain.
                  //Thus for the first timestep we calculate an intial viscosity based on the strain rate. Why not later: seems to cause unstable solution(?)
                  if (this->get_timestep_number() == 1)
                    {
                      const double edot_ii=std::max(std::sqrt(std::max(-second_invariant(deviator(strain_rate)), 0.)),
                                                    min_strain_rate);
                      scalar_viscosity= 1/Gamma *  std::pow(edot_ii,((1. - n)/n));
                    }

                  double n_iterations = 1;
                  double max_iteration = 100;
                  double residual = scalar_viscosity;
                  double threshold = 0.0001*scalar_viscosity;
                  SymmetricTensor<2,dim> stress;
                  stress = scalar_viscosity * V_r4 * deviatoric_strain_rate / 1e6; // Use stress in MPa
                  // std::cout << "Initial stress: " << stress << std::endl;
                  while (std::abs(residual) > threshold && n_iterations < max_iteration)
                    // while (n_iterations < max_iteration)
                    {
                      // std::cout << "n_iterations: " << n_iterations << std::endl;
                      stress = (1./2.) * (stress + scalar_viscosity * V_r4 * deviatoric_strain_rate / 1e6);
                      // std::cout << "old_stress_strain_director " << old_stress_strain_director << std::endl;
                      // std::cout << "deviatoric_strain_rate " << deviatoric_strain_rate << std::endl;
                      // std::cout << "Anisotropic stress " << stress << std::endl;

                      Tensor<2,3> S_CPO=transpose(R)*stress*R;
                      // std::cout << "stress " << stress <<std::endl;
                      // std::cout << "stress CPO " << S_CPO <<std::endl;

                      double Jhill = F*pow((S_CPO[0][0]-S_CPO[1][1]),2) + G*pow((S_CPO[1][1]-S_CPO[2][2]),2) + H*pow((S_CPO[2][2]-S_CPO[0][0]),2) + 2*L*pow(S_CPO[1][2],2) + 2*M*pow(S_CPO[0][2],2) + 2*N*pow(S_CPO[0][1],2);
                      if (Jhill < 0)
                        {
                          Jhill = std::abs(F)*pow((S_CPO[0][0]-S_CPO[1][1]),2) + std::abs(G)*pow((S_CPO[1][1]-S_CPO[2][2]),2) + std::abs(H)*pow((S_CPO[2][2]-S_CPO[0][0]),2) + 2*L*pow(S_CPO[1][2],2) + 2*M*pow(S_CPO[0][2],2) + 2*N*pow(S_CPO[0][1],2);
                        }
                      // std::cout << "Jhill " << Jhill <<std::endl;

                      AssertThrow(isfinite(Jhill),
                                  ExcMessage("Jhill should be finite"));
                      AssertThrow(Jhill >= 0,
                                  ExcMessage("Jhill should not be negative"));

                      double scalar_viscosity_new = (1 / (Gamma * std::pow(Jhill,(n-1)/2))) * 1e6; // convert from MPa to Pa
                      residual = std::abs(scalar_viscosity_new - scalar_viscosity);
                      scalar_viscosity = scalar_viscosity_new;
                      threshold = 0.001*scalar_viscosity;
                      // std::cout << "scalar_viscosity in loop: " << scalar_viscosity <<std::endl;
                      // std::cout << "residual: " << residual <<std::endl;
                      n_iterations += 1;

                      // std::cout<<"F "<<F<<" G "<<G<<" H "<<H<<" L "<<L<<" M "<<M<<" N "<<N<<std::endl;
                      // std::cout << "old_stress_strain_director " << old_stress_strain_director << std::endl;
                      // std::cout << "deviatoric_strain_rate " << deviatoric_strain_rate << std::endl;
                      // std::cout << "stress " << stress <<std::endl;
                      // std::cout << "R " << R <<std::endl;
                      // std::cout << "stress CPO " << S_CPO <<std::endl;
                      // std::cout << "Jhill " << Jhill <<std::endl;
                    }
                  //Overwrite the scalar viscosity with an effective viscosity
                  out.viscosities[q] = scalar_viscosity;//composition[ind_vis];//
                  //std::cout << "Final scalar_viscosity: " << scalar_viscosity <<std::endl;

                  AssertThrow(out.viscosities[q] > 0,
                              ExcMessage("Viscosity should be positive"));
                  AssertThrow(isfinite(out.viscosities[q]),
                              ExcMessage("Viscosity should be finite"));
                  //Compute Rotation matrix

                }
            }
          else
            {
              // std::cout << "AV is nullptr since dev_strainrate is 0: " << determinant(deviatoric_strain_rate) << std::endl;
              if (anisotropic_viscosity != nullptr)
                {
                  // std::cout << "AV is not nullptr but isotropic " << std::endl;
                  SymmetricTensor<2,6> V;
                  V[0][0] = 2.0/3.0;
                  V[0][1] = -1.0/3.0;
                  V[0][2] = -1.0/3.0;
                  V[1][1] = 2.0/3.0;
                  V[1][2] = -1.0/3.0;
                  V[2][2] = 2.0/3.0;
                  V[3][3] = 1;
                  V[4][4] = 1;
                  V[5][5] = 1;
                  //Convert rank 2 viscosity tensor to rank 4
                  FullMatrix<double> V_mat(6,6);
                  for (unsigned int vi=0; vi<6; ++vi)
                    {
                      for (unsigned int vj=0; vj<6; ++vj)
                        {
                          V_mat[vi][vj] = V[vi][vj];
                        }
                    }
                  SymmetricTensor<4,dim> V_r4;
                  dealii::Physics::Notation::Kelvin::to_tensor(V_mat, V_r4);
                  anisotropic_viscosity->stress_strain_directors[q] = V_r4;

                }
            }
          // Prescribe the stress strain directors and scalar viscosity to compositional field for access in the next time step
          if (PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>())
            {
              //std::vector<double> ViscoTensor_array(SymmetricTensor<4,dim>::n_independent_components);
              // FullMatrix<double> V_mat = dealii::Physics::Notation::Kelvin::to_matrix(anisotropic_viscosity->stress_strain_directors[q]);
              // SymmetricTensor<2,6> V_r2;
              // for (unsigned int vi=0; vi<6; ++vi)
              //   {
              //     for (unsigned int vj=0; vj<6; ++vj)
              //       {
              //         V_r2[vi][vj] = V_mat[vi][vj];
              //       }
              //   }
              //std::copy(anisotropic_viscosity->stress_strain_directors[q].begin_raw(), anisotropic_viscosity->stress_strain_directors[q].end_raw(), ViscoTensor_array.begin());
              //for (unsigned int i = 0; i < SymmetricTensor<4,dim>::n_independent_components ; ++i)
              //{
              //const unsigned int ind = this->introspection().compositional_index_for_name(ssd_names[i]);
              //prescribed_field_out->prescribed_field_outputs[q][ind] = ViscoTensor_array[i];
              //AssertThrow(isfinite(ViscoTensor_array[i]),
              //ExcMessage("Assigning prescribed field should be finite"));
              //}
              const unsigned int ind_vis = this->introspection().compositional_index_for_name("scalar_vis");
              prescribed_field_out->prescribed_field_outputs[q][ind_vis] = out.viscosities[q];
              // std::cout << "Saved ViscoTensor_array: ";
              // for (double VT_a_item: ViscoTensor_array)
              //     std::cout << VT_a_item << ' ';
              // std::cout << std::endl;
              // std::cout << "Saved ssd: " << anisotropic_viscosity->stress_strain_directors[q] << std::endl;
              //std::cout << "Saved viscosity: " << out.viscosities[q] << std::endl;
            }
        }
    }



    template <int dim>
    bool
    LPO_AV_3D<dim>::is_compressible () const
    {
      return false;
    }



    // template <int dim>
    // double
    // LPO_AV_3D<dim>::reference_viscosity () const
    // {
    //   return 1e20;
    // }



    template <int dim>
    void
    LPO_AV_3D<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AV Hill");
        {

          equation_of_state.parse_parameters (prm);
          eta = prm.get_double("Reference viscosity");
          min_strain_rate = prm.get_double("Minimum strain rate");
          grain_size = prm.get_double("Grain size");
          CnI_F = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for F")));
          CnI_G = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for G")));
          CnI_H = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for H")));
          CnI_L = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for L")));
          CnI_M = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for M")));
          CnI_N = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(prm.get("Coefficients and intercept for N")));
          compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // prm.enter_subsection("Particles");
      // {
      //   prm.enter_subsection("Crystal Preferred Orientation");
      //   {
      //     n_grains = prm.get_integer("Number of grains per particle");
      //   }
      //   prm.leave_subsection();
      // }
      // prm.leave_subsection();

      // Declare dependence
      this->model_dependence.density = NonlinearDependence::compositional_fields;
    }



    template <int dim>
    void
    LPO_AV_3D<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AV Hill");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm);
          // prm.declare_entry ("Coefficients and intercept for F", "4.600447869251792, -3.432107596, -0.132855296, -0.00130425,  2.360195770340767,  -1.503338915, 1.2635426953453306, 0.0071666,  1.2241252112687866, -2.655458732, -1.734603653, 0.039704515, 2.241765907780337",
          //                    Patterns::List(Patterns::Double()),
          //                    "6 Coefficients and 1 intercept to compute the Hill Parameter F.");
          // prm.declare_entry ("Coefficients and intercept for G", "-4.684144376,  2.4501679747908645, -1.291945979, -0.000447184, 0.7364410668996572, -1.134353994, -1.369971542, 0.006638505,  1.8100476263969771, 0.278455287,  2.0422841708796087, -0.016914332, 0.37853585511038185",
          //                    Patterns::List(Patterns::Double()),
          //                    "6 Coefficients and 1 intercept to compute the Hill Parameter G.");
          // prm.declare_entry ("Coefficients and intercept for H", "4.020574320458808, -2.637132566, 0.5755508961574253, 0.002686342,  1.5151533596671836, 0.3489226631973101, 1.3802207762705025, -0.007362419, 1.4996869664711834, -4.044348002, -1.978212652, 0.056659336732010734, 0.37853585511038185",
          //                    Patterns::List(Patterns::Double()),
          //                    "6 Coefficients and 1 intercept to compute the Hill Parameter H.");
          // prm.declare_entry ("Coefficients and intercept for L", "1.3973326672504107,  -1.81296565,  0.8981686897139447, -0.001066616, -0.063667124, -0.374885301, 0.7288437733735754, -0.000503799, -0.99360894,  -0.153230898, -1.019483357, 0.007165858, 1.9767970726829662",
          //                    Patterns::List(Patterns::Double()),
          //                    "6 Coefficients and 1 intercept to compute the Hill Parameter L.");
          // prm.declare_entry ("Coefficients and intercept for M", "1.6321611239635443,  0.30732770773835577,  1.1935547119876715, 0.002980557,  3.0859826232991554, -0.954888925, 0.6077177400767578, -0.009972845, -3.089290395, 0.648574284,  -0.558830801, -0.007724845, 0.8743203613069472",
          //                    Patterns::List(Patterns::Double()),
          //                    "6 Coefficients and 1 intercept to compute the Hill Parameter M.");
          // prm.declare_entry ("Coefficients and intercept for N", "0.558022129, 0.27369095729305215,  0.11267428759159648,  0.001899995,  -1.59528053,  0.17253325907449615,  0.4514508967688472, -0.010286343, 2.9328835681486245, -1.375833649, -0.045106095, 0.015508168775970624, 1.338955424846116",
          //                    Patterns::List(Patterns::Double()),
          //                    "6 Coefficients and 1 intercept to compute the Hill Parameter N.");

          // old
          prm.declare_entry ("Coefficients and intercept for F", "1.0390459583037057,  -0.767458622,  0.003066208,  0.19651133418307049,  0.413093763,  0.015463162,  -0.935925291,  -2.392877563,  0.051834768, 1.0799807050187482",
                             Patterns::List(Patterns::Double()),
                             "6 Coefficients and 1 intercept to compute the Hill Parameter F.");
          prm.declare_entry ("Coefficients and intercept for G", "-2.836270315,  -1.632453092,  0.000687606,  0.2671850239576621,  -0.993392913,  0.002699241,  1.9689530759060374,  2.314442451425019,  -0.018655905, 0.6887411607403755",
                             Patterns::List(Patterns::Double()),
                             "6 Coefficients and 1 intercept to compute the Hill Parameter G.");
          prm.declare_entry ("Coefficients and intercept for H", "1.6687493021559732,  0.5797579293682223,  0.003241593,  0.701661336,  0.2513824481429968,  0.000229291,  -2.003227619,  -2.57032429,  0.071454541, 0.7490268673620638",
                             Patterns::List(Patterns::Double()),
                             "6 Coefficients and 1 intercept to compute the Hill Parameter H.");
          prm.declare_entry ("Coefficients and intercept for L", "-0.325145943,  0.7284642859944138,  0.000404879,  -0.665446098,  0.5152847961409479,  0.002722782,  -1.026786493,  -1.262574542,  0.009168498, 1.595422603",
                             Patterns::List(Patterns::Double()),
                             "6 Coefficients and 1 intercept to compute the Hill Parameter L.");
          prm.declare_entry ("Coefficients and intercept for M", "1.6427437063774875,  0.8777500120437522,  0.004651732,  2.489417876177839,  0.8162729707609052,  -0.010736521,  -2.49420455,  -0.511446494,  -0.009362491, 0.893677343",
                             Patterns::List(Patterns::Double()),
                             "6 Coefficients and 1 intercept to compute the Hill Parameter M.");
          prm.declare_entry ("Coefficients and intercept for N", "0.8122098589701904,  0.15663795996228266,  0.001500252,  -1.648578168,  0.19362392490527092,  -0.009650519,  1.6796559729985163,  -0.103640482,  0.01971017, 1.2132200780065174",
                             Patterns::List(Patterns::Double()),
                             "6 Coefficients and 1 intercept to compute the Hill Parameter N.");

          prm.declare_entry ("Reference viscosity", "1e9",
                             Patterns::Double(),
                             "Magnitude of reference viscosity.");
          prm.declare_entry ("Minimum strain rate", "1.4e-20", Patterns::Double(),
                             "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}");
          prm.declare_entry ("Grain size", "1000",
                             Patterns::Double(),
                             "Olivine anisotropic viscosity is dependent of grain size. Value is given in microns");
          prm.declare_entry ("Density differential for compositional field 1", "0.",
                             Patterns::Double(),
                             "If compositional fields are used, then one would frequently want "
                             "to make the density depend on these fields. In this simple material "
                             "model, we make the following assumptions: if no compositional fields "
                             "are used in the current simulation, then the density is simply the usual "
                             "one with its linear dependence on the temperature. If there are compositional "
                             "fields, then the density only depends on the first one in such a way that "
                             "the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
                             "This parameter describes the value of $\\Delta \\rho$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}/unit change in composition.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }



    template <int dim>
    void
    LPO_AV_3D<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<AV<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::AV<dim>> (n_points));
        }

      if (out.template get_additional_output<PrescribedFieldOutputs<dim>>() == NULL)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points,this->n_compositional_fields()));
        }
    }
  }
}



// explicit instantiations
namespace aspect
{
  namespace Assemblers
  {
#define INSTANTIATE(dim) \
  template class StokesPreconditionerAV<dim>; \
  template class StokesIncompressibleTermsAV<dim>; \
  //template class StokesBoundaryTractionAV<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingAV,
                                  "anisotropic shear heating for LPO_AV_3D",
                                  "Implementation of a standard model for shear heating. "
                                  "Adds the term: "
                                  "$  2 \\eta \\left( \\varepsilon - \\frac{1}{3} \\text{tr} "
                                  "\\varepsilon \\mathbf 1 \\right) : \\left( \\varepsilon - \\frac{1}{3} "
                                  "\\text{tr} \\varepsilon \\mathbf 1 \\right)$ to the "
                                  "right-hand side of the temperature equation.")
  }



  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(LPO_AV_3D,
                                   "LPO Anisotropic Viscosity Hill material",
                                   "Olivine LPO related viscous anisotropy based on the Simple material model")
  }
}
}
