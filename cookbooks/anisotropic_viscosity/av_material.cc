/*
 Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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

#include <aspect/introspection.h>
#include <aspect/material_model/interface.h>
#include <aspect/plugins.h>
#include <aspect/simulator/assemblers/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#ifndef __aspect__av_material_h
#define __aspect__av_material_h

#include <aspect/simulator_access.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/base/geometry_info.h>
#include <aspect/simulator_access.h>

#include <aspect/material_model/simple.h>
#include <aspect/material_model/grain_size.h>
#include <aspect/heating_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator_signals.h>
#include <aspect/postprocess/particles.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>

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
    template <int dim>
    class AnisotropicViscosity : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        AnisotropicViscosity(const unsigned int n_points);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Stress-strain "director" tensors at the given positions. This
         * variable is used to implement anisotropic viscosity.
         *
         * @note The strain rate term in equation (1) of the manual will be
         * multiplied by this tensor *and* the viscosity scalar ($\eta$), as
         * described in the manual section titled "Constitutive laws". This
         * variable is assigned the rank-four identity tensor by default.
         * This leaves the isotropic constitutive law unchanged if the material
         * model does not explicitly assign a value.
         */
        std::vector<SymmetricTensor<4,dim> > stress_strain_directors;
    };

    namespace
    {



      template <int dim>
      std::vector<std::string> make_AnisotropicViscosity_additional_outputs_names()
      {
        std::vector<std::string> names;

        for (unsigned int i = 0; i < Tensor<4,dim>::n_independent_components ; ++i)
          {
            TableIndices<4> indices(Tensor<4,dim>::unrolled_to_component_indices(i));
            names.push_back("anisotropic_viscosity"+std::to_string(indices[0])+std::to_string(indices[1])+std::to_string(indices[2])+std::to_string(indices[3]));
          }
        return names;
      }
    }



    template <int dim>
    AnisotropicViscosity<dim>::AnisotropicViscosity (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(make_AnisotropicViscosity_additional_outputs_names<dim>()),
      stress_strain_directors(n_points, dealii::identity_tensor<dim> ())
    {}



    template <int dim>
    std::vector<double>
    AnisotropicViscosity<dim>::get_nth_output(const unsigned int idx) const
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
    class StokesPreconditionerAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        /**
         * Create AnisotropicViscosities.
         */
        virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };

    /**
     * This class assembles the terms for the matrix and right-hand-side of the incompressible
     * Stokes equation for the current cell.
     */
    template <int dim>
    class StokesIncompressibleTermsAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;

        /**
         * Create AdditionalMaterialOutputsStokesRHS if we need to do so.
         */
        virtual void create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;
    };



    template <int dim>
    void
    StokesPreconditionerAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesPreconditioner<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesPreconditioner<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesPreconditioner<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesPreconditioner<dim>& > (data_base);

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

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
    StokesPreconditionerAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }
    }



    template <int dim>
    void
    StokesIncompressibleTermsAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
      const double pressure_scaling = this->get_pressure_scaling();

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = scratch.material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >();

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
    StokesIncompressibleTermsAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      const unsigned int n_points = outputs.viscosities.size();

      if (outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      if (this->get_parameters().enable_additional_stokes_rhs
          && outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >() == nullptr)
        {
          outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));
        }
      Assert(!this->get_parameters().enable_additional_stokes_rhs
             ||
             outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
             == n_points, ExcInternalError());
    }
  }

  namespace HeatingModel
  {
    template <int dim>
    class ShearHeatingAnisotropicViscosity : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * Allow the heating model to attach additional material model outputs.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const;
    };



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.strain_rate.size(),
             ExcMessage ("The shear heating plugin needs the strain rate!"));

      // Some material models provide dislocation viscosities and boundary area work fractions
      // as additional material outputs. If they are attached, use them.
      const MaterialModel::DislocationViscosityOutputs<dim> *disl_viscosities_out =
        material_model_outputs.template get_additional_output<MaterialModel::DislocationViscosityOutputs<dim> >();

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

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

          // If dislocation viscosities and boundary area work fractions are provided, reduce the
          // overall heating by this amount (which is assumed to increase surface energy)
          if (disl_viscosities_out != 0)
            {
              heating_model_outputs.heating_source_terms[q] *= 1 - disl_viscosities_out->boundary_area_change_work_fractions[q] *
                                                               material_model_outputs.viscosities[q] / disl_viscosities_out->dislocation_viscosities[q];
            }

          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }



    template <int dim>
    void
    ShearHeatingAnisotropicViscosity<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      const unsigned int n_points = material_model_outputs.viscosities.size();

      if (material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >() == nullptr)
        {
          material_model_outputs.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
        }

      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }
  }

  namespace MaterialModel
  {
    // The AV material model calculates an anisotropic viscosity tensor from director vectors and the normal and shear
    // viscosities (defined in the .prm file). In contrast to the `Anisotropic` material model, the principal directions of the
    // tensor are not read from an input file, but instead are computed at every quadrature point. This material model is used in T independent models.
    // All the parameters are defined below.
    template <int dim>
    class AV : public MaterialModel::Simple<dim>
    {
      public:
        virtual void initialize();
        virtual void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const;
        static void declare_parameters (ParameterHandler &prm);
        virtual void parse_parameters (ParameterHandler &prm);
        virtual bool is_compressible () const;
        virtual double reference_viscosity () const;
        virtual double reference_density () const;
        virtual void create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const;
      private:
        double eta_N, viscosity_ratio; //normal viscosity and ratio between the shear and the normal viscosities
        static int delta (const unsigned int i, const unsigned int j); //kronecker delta function
        void set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const;
    };
  }
}

namespace aspect
{

//Next session is a more evolved implementation of anisotropic viscosity in the material model based on Jonathan Perry-Houts'paper
  namespace MaterialModel
  {
    template <int dim>
    void
    AV<dim>::set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const
    {
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesPreconditioner<dim>>(*(assemblers.stokes_preconditioner[i])))
            assemblers.stokes_preconditioner[i] = std_cxx14::make_unique<Assemblers::StokesPreconditionerAnisotropicViscosity<dim> > ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (Plugins::plugin_type_matches<Assemblers::StokesIncompressibleTerms<dim>>(*(assemblers.stokes_system[i])))
            assemblers.stokes_system[i] = std_cxx14::make_unique<Assemblers::StokesIncompressibleTermsAnisotropicViscosity<dim> > ();
        }
    }



    template <int dim>
    void
    AV<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&AV<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
      AssertThrow((dim==2),
                  ExcMessage("For now the anisotropic viscosity assemblers work only in 2D"));

    }



    template <int dim>
    void
    AV<dim>::evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        out.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      AssertThrow((this->introspection().compositional_name_exists("gamma")),
                  ExcMessage("AV material model only works if there is a compositional field called gamma."));
      AssertThrow(this->introspection().compositional_name_exists("ni"),
                  ExcMessage("AV material model only works if there is a compositional field called ni."));
      AssertThrow(this->introspection().compositional_name_exists("nj"),
                  ExcMessage("AV material model only works if there is a compositional field called nj."));
      if (dim == 3)
        AssertThrow(this->introspection().compositional_name_exists("nk"),
                    ExcMessage("AV material model only works if there is a compositional field called nk."));

      const unsigned int c_idx_gamma = this->introspection().compositional_index_for_name("gamma");

      std::vector<unsigned int> c_idx_n;
      c_idx_n.push_back (this->introspection().compositional_index_for_name("ni"));
      c_idx_n.push_back (this->introspection().compositional_index_for_name("nj"));
      if (dim == 3)
        c_idx_n.push_back (this->introspection().compositional_index_for_name("nk"));

      // Get the grad_u tensor, at the center of this cell, if possible.

      std::vector<Tensor<2,dim> > velocity_gradients (in.n_evaluation_points());
      if (in.current_cell.state() == IteratorState::valid)
        {
          std::vector<Point<dim> > quadrature_positions(in.n_evaluation_points());
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            quadrature_positions[i] = this->get_mapping().transform_real_to_unit_cell(in.current_cell, in.position[i]);

          // FEValues requires a quadrature and we provide the default quadrature
          // as we only need to evaluate the gradients of the solution.
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   Quadrature<dim>(quadrature_positions),
                                   update_gradients);
          fe_values.reinit (in.current_cell);
          fe_values[this->introspection().extractors.velocities]
          .get_function_gradients(this->get_solution(), velocity_gradients);
        }

      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          out.densities[q] = (in.composition[q][c_idx_gamma] > 0.8 ? 1 : 0);
          out.viscosities[q] = eta_N;
          out.thermal_expansion_coefficients[q] = 0;
          out.specific_heat[q] = 0;
          out.thermal_conductivities[q] = 0;
          out.compressibilities[q] = 0.0;
          out.entropy_derivative_pressure[q] = 0.0;
          out.entropy_derivative_temperature[q] = 0.0;
          for (unsigned int c=0; c<in.composition[q].size(); ++c)
            out.reaction_terms[q][c] = 0.0;

          Tensor<1,dim> n;
          for (unsigned int i=0; i<dim; ++i)
            n[i] = in.composition[q][c_idx_n[i]];

          // The computation of the viscosity tensor is only
          // necessary after the simulator has been initialized.
          if (this->simulator_is_past_initialization())
            {
              Tensor<1,dim> n_dot;
              if (n.norm() > 0.5 && in.composition[q][c_idx_gamma] > 0.8)
                {
                  // Symmetric and anti-symmetric parts of grad_u
                  const SymmetricTensor<2,dim> D = symmetrize(velocity_gradients[q]);
                  Tensor<2,dim> W;
                  for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      W[i][j] = velocity_gradients[q][i][j] - D[i][j];

                  for (unsigned int i=0; i<dim; ++i)
                    {
                      // outer summation over each value of j.
                      for (unsigned int j=0; j<dim; ++j)
                        {
                          float Wn = W[i][j];
                          for (unsigned int k=0; k<dim; ++k)
                            {
                              Wn -= D[k][i]*n[k]*n[j] - D[k][j]*n[k]*n[i];
                            }
                          n_dot[i] += Wn * n[j];
                        }
                    }

                  // make sure that n is a unit vector. We have checked
                  // above that n != 0.0.
                  // handle time step 0 differently, because time step length is 0
                  if (this->get_timestep() == 0)
                    {
                      n_dot = (n/n.norm())-n;
                    }
                  else
                    {
                      Tensor<1,dim> n_new = n + (n_dot * this->get_timestep());
                      Assert (n_new.norm() != 0, ExcInternalError());
                      n_new /= n_new.norm();
                      n_dot = (n_new-n)/ this->get_timestep();
                    }
                }
              else
                {
                  n_dot = 0;
                }
              // update  n[i] = in.composition[q][c_idx_n[i]] with adding to it n_dot*dt
              for (unsigned int i=0; i<dim; ++i)
                if (this->get_timestep() == 0)
                  out.reaction_terms[q][c_idx_n[i]] = n_dot[i];
                else
                  out.reaction_terms[q][c_idx_n[i]] = n_dot[i] * this->get_timestep();

              if (n.norm() > 0.5)
                {
                  n /= n.norm();
                  SymmetricTensor<4,dim> Lambda;
                  for (unsigned int i=0; i<dim; ++i)
                    for (unsigned int j=0; j<dim; ++j)
                      for (unsigned int k=0; k<dim; ++k)
                        for (unsigned int l=0; l<dim; ++l)
                          Lambda[i][j][k][l] = 1./2. * (n[i]*n[k]*delta(l,j)
                                                        + n[j]*n[k]*delta(i,l)
                                                        + n[i]*n[l]*delta(k,j)
                                                        + n[j]*n[l]*delta(i,k))
                                               - 2*n[i]*n[j]*n[k]*n[l];

                  if (anisotropic_viscosity != nullptr)
                    {
                      anisotropic_viscosity->stress_strain_directors[q] =dealii::identity_tensor<dim> ()
                                                                         - (1. - viscosity_ratio) * Lambda;
                      SymmetricTensor<2,3> ViscoTensor;
                      ViscoTensor[0][0]=anisotropic_viscosity->stress_strain_directors[q][0][0][0][0];
                      ViscoTensor[0][1]=anisotropic_viscosity->stress_strain_directors[q][0][0][1][1];
                      ViscoTensor[0][2]=anisotropic_viscosity->stress_strain_directors[q][0][0][0][1] * std::sqrt(2);
                      ViscoTensor[1][0]=anisotropic_viscosity->stress_strain_directors[q][1][1][0][0];
                      ViscoTensor[1][1]=anisotropic_viscosity->stress_strain_directors[q][1][1][1][1];
                      ViscoTensor[1][2]=anisotropic_viscosity->stress_strain_directors[q][1][1][0][1] * std::sqrt(2);
                      ViscoTensor[2][0]=anisotropic_viscosity->stress_strain_directors[q][0][1][0][0] * std::sqrt(2);
                      ViscoTensor[2][1]=anisotropic_viscosity->stress_strain_directors[q][0][1][1][1] * std::sqrt(2);
                      ViscoTensor[2][2]=anisotropic_viscosity->stress_strain_directors[q][0][1][0][1] * 2;

                      const std::array<double,3> Viscoeigenvalues = eigenvalues(ViscoTensor);
                      for (unsigned int i=0; i<3; ++i)
                        {
                          AssertThrow((Viscoeigenvalues[i]>0),
                                      ExcMessage("Eigenvalue "+ std::to_string(i) + " of the viscosity tensor is negative at " +
                                                 std::to_string(Viscoeigenvalues[i]) + ". This is not allowed."));
                        }
                    }
                }
            }
        }
    }



    template <int dim>
    int
    AV<dim>::delta (const unsigned int i,
                    const unsigned int j)
    {
      return (i == j ? 1 : 0);
    }



    template <int dim>
    bool
    AV<dim>::is_compressible () const
    {
      return false;
    }



    template <int dim>
    double
    AV<dim>::reference_density () const
    {
      return 1.0;
    }



    template <int dim>
    double
    AV<dim>::reference_viscosity () const
    {
      return 1.0;
    }



    template <int dim>
    void
    AV<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AV");
        {
          eta_N = prm.get_double("Normal viscosity");
          viscosity_ratio = prm.get_double("Shear viscosity")/eta_N;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AV<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("AV");
        {
          prm.declare_entry ("Normal viscosity", "1e-3",
                             Patterns::Double(),
                             "Magnitude of normal viscosity.");
          prm.declare_entry ("Shear viscosity", "1e-4",
                             Patterns::Double(),
                             "Magnitude of shear viscosity.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AV<dim>::create_additional_named_outputs(MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<AnisotropicViscosity<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AnisotropicViscosity<dim>> (n_points));
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
  template class StokesPreconditioner<dim>; \
  template class StokesIncompressibleTerms<dim>; \
  template class StokesBoundaryTraction<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ShearHeatingAnisotropicViscosity,
                                  "anisotropic shear heating",
                                  "Implementation of a standard model for shear heating. "
                                  "Adds the term: "
                                  "$  2 \\eta \\left( \\varepsilon - \\frac{1}{3} \\text{tr} "
                                  "\\varepsilon \\mathbf 1 \\right) : \\left( \\varepsilon - \\frac{1}{3} "
                                  "\\text{tr} \\varepsilon \\mathbf 1 \\right)$ to the "
                                  "right-hand side of the temperature equation.")
  }



  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AV,
                                   "AV material",
                                   "Transverse isotropic material model.")
  }
}
#endif
