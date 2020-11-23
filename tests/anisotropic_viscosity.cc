#include <aspect/simulator_access.h>

#include <aspect/material_model/simple.h>
#include <aspect/material_model/grain_size.h>
#include <aspect/heating_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/simulator_signals.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace MaterialModel
  {
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

    template <int dim>
    AnisotropicViscosity<dim>::AnisotropicViscosity (const unsigned int n_points)
      :
      NamedAdditionalMaterialOutputs<dim>(std::vector<std::string>(1,"anisotropic_viscosity")),
      stress_strain_directors(n_points, dealii::identity_tensor<dim> ())
    {}



    template <int dim>
    std::vector<double>
    AnisotropicViscosity<dim>::get_nth_output(const unsigned int /*idx*/) const
    {
      return std::vector<double>();
    }
  }

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
     * A class containing the functions to assemble the compressible adjustment
     * to the Stokes preconditioner.
     */
    template <int dim>
    class StokesCompressiblePreconditionerAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
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

    /**
     * This class assembles the term that arises in the viscosity term of Stokes matrix for
     * compressible models, because the divergence of the velocity is not longer zero.
     */
    template <int dim>
    class StokesCompressibleStrainRateViscosityTermAnisotropicViscosity : public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        execute(internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                internal::Assembly::CopyData::CopyDataBase<dim> &data) const;
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

          const bool use_tensor = (anisotropic_viscosity != nullptr);

          const SymmetricTensor<4, dim> &stress_strain_director = (use_tensor)
                                                                  ?
                                                                  anisotropic_viscosity->stress_strain_directors[q]
                                                                  :
                                                                  dealii::identity_tensor<dim>();



          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += ((
                                              use_tensor ?
                                              2.0 * eta * (scratch.grads_phi_u[i]
                                                           * stress_strain_director
                                                           * scratch.grads_phi_u[j]) :
                                              2.0 * eta * (scratch.grads_phi_u[i]
                                                           * scratch.grads_phi_u[j]))
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
    StokesCompressiblePreconditionerAnisotropicViscosity<dim>::
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
                  scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);

                  ++i_stokes;
                }
              ++i;
            }

          const double eta_two_thirds = scratch.material_model_outputs.viscosities[q] * 2.0 / 3.0;

          const bool use_tensor = (anisotropic_viscosity != nullptr);

          const SymmetricTensor<4, dim> &stress_strain_director = (use_tensor)
                                                                  ?
                                                                  anisotropic_viscosity->stress_strain_directors[q]
                                                                  :
                                                                  dealii::identity_tensor<dim>();

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i = 0; i < stokes_dofs_per_cell; ++i)
            for (unsigned int j = 0; j < stokes_dofs_per_cell; ++j)
              if (scratch.dof_component_indices[i] ==
                  scratch.dof_component_indices[j])
                data.local_matrix(i, j) += (- (use_tensor ?
                                               eta_two_thirds * (scratch.div_phi_u[i] * trace(stress_strain_director * scratch.grads_phi_u[j]))
                                               :
                                               eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                              ))
                                           * JxW;
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

          const bool use_tensor = (anisotropic_viscosity != nullptr);

          const SymmetricTensor<4, dim> &stress_strain_director = (use_tensor)
                                                                  ?
                                                                  anisotropic_viscosity->stress_strain_directors[q]
                                                                  :
                                                                  dealii::identity_tensor<dim>();

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
                    data.local_matrix(i,j) += ( (use_tensor ?
                                                 eta * 2.0 * (scratch.grads_phi_u[i] * stress_strain_director * scratch.grads_phi_u[j])
                                                 :
                                                 eta * 2.0 * (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]))
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



    template <int dim>
    void
    StokesCompressibleStrainRateViscosityTermAnisotropicViscosity<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      if (!scratch.rebuild_stokes_matrix)
        return;

      const MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        scratch.material_model_outputs.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = this->get_fe();
      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
            {
              if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                {
                  scratch.grads_phi_u[i_stokes] = scratch.finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                  scratch.div_phi_u[i_stokes]   = scratch.finite_element_values[introspection.extractors.velocities].divergence (i, q);

                  ++i_stokes;
                }
              ++i;
            }

          // Viscosity scalar
          const double eta_two_thirds = scratch.material_model_outputs.viscosities[q] * 2.0 / 3.0;

          const bool use_tensor = (anisotropic_viscosity != nullptr);

          const SymmetricTensor<4, dim> &stress_strain_director = (use_tensor)
                                                                  ?
                                                                  anisotropic_viscosity->stress_strain_directors[q]
                                                                  :
                                                                  dealii::identity_tensor<dim>();

          const double JxW = scratch.finite_element_values.JxW(q);

          for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
            for (unsigned int j=0; j<stokes_dofs_per_cell; ++j)
              {
                data.local_matrix(i,j) += (- (use_tensor ?
                                              eta_two_thirds * (scratch.div_phi_u[i] * trace(stress_strain_director * scratch.grads_phi_u[j]))
                                              :
                                              eta_two_thirds * (scratch.div_phi_u[i] * scratch.div_phi_u[j])
                                             ))
                                          * JxW;
              }
        }
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
    using namespace dealii;

    template <int dim>
    class Anisotropic : public MaterialModel::Simple<dim>
    {
      public:
        virtual void initialize();

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
        static void declare_parameters(ParameterHandler &prm);
        virtual void parse_parameters(ParameterHandler &prm);

        /**
          * Return true if the compressibility() function returns something that
          * is not zero.
          */
        virtual bool
        is_compressible () const;

      private:
        void set_assemblers(const SimulatorAccess<dim> &,
                            Assemblers::Manager<dim> &assemblers) const;

        // Constitutive tensor
        SymmetricTensor<4,dim> C;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Anisotropic<dim>::set_assemblers(const SimulatorAccess<dim> &,
                                     Assemblers::Manager<dim> &assemblers) const
    {
      for (unsigned int i=0; i<assemblers.stokes_preconditioner.size(); ++i)
        {
          if (dynamic_cast<Assemblers::StokesPreconditioner<dim> *>(assemblers.stokes_preconditioner[i].get()) != nullptr)
            assemblers.stokes_preconditioner[i] = std_cxx14::make_unique<Assemblers::StokesPreconditionerAnisotropicViscosity<dim> > ();
          if (dynamic_cast<Assemblers::StokesCompressiblePreconditioner<dim> *>(assemblers.stokes_preconditioner[i].get()) != nullptr)
            assemblers.stokes_preconditioner[i] = std_cxx14::make_unique<Assemblers::StokesCompressiblePreconditionerAnisotropicViscosity<dim> > ();
        }

      for (unsigned int i=0; i<assemblers.stokes_system.size(); ++i)
        {
          if (dynamic_cast<Assemblers::StokesIncompressibleTerms<dim> *>(assemblers.stokes_system[i].get()) != nullptr)
            assemblers.stokes_system[i] = std_cxx14::make_unique<Assemblers::StokesIncompressibleTermsAnisotropicViscosity<dim> > ();
          if (dynamic_cast<Assemblers::StokesCompressibleStrainRateViscosityTerm<dim> *>(assemblers.stokes_system[i].get()) != nullptr)
            assemblers.stokes_system[i] = std_cxx14::make_unique<Assemblers::StokesCompressibleStrainRateViscosityTermAnisotropicViscosity<dim> > ();
        }
    }

    template <int dim>
    void
    Anisotropic<dim>::
    initialize()
    {
      this->get_signals().set_assemblers.connect (std::bind(&Anisotropic<dim>::set_assemblers,
                                                            std::cref(*this),
                                                            std::placeholders::_1,
                                                            std::placeholders::_2));
    }

    template <int dim>
    void
    Anisotropic<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      MaterialModel::AnisotropicViscosity<dim> *anisotropic_viscosity =
        out.template get_additional_output<MaterialModel::AnisotropicViscosity<dim> >();

      Simple<dim>::evaluate(in, out);
      Point<dim> center;
      center[0] = 0.5;
      center[1] = 0.5;
      if (dim == 3)
        center[2] = 0.5;
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          const double pressure = in.pressure[i];
          out.densities[i] = 1.0 + pressure;
          out.compressibilities[i] = 1.0 / (1. + pressure);
          if ((in.position[i]-center).norm() < 0.25)
            {
              if (anisotropic_viscosity != nullptr)
                anisotropic_viscosity->stress_strain_directors[i] = C;
            }
        }
    }

    template <int dim>
    bool
    Anisotropic<dim>::
    is_compressible () const
    {
      return true;
    }

    template <int dim>
    void
    Anisotropic<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      Simple<dim>::declare_parameters (prm);
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Anisotropic");
        {
          if (dim == 2)
            prm.declare_entry ("Viscosity tensor",
                               "1, 0, 0,"
                               "0, 1, 0,"
                               "0, 0,.5",
                               Patterns::List(Patterns::Double()),
                               "Viscosity-scaling tensor in Voigt notation.");
          else
            prm.declare_entry ("Viscosity tensor",
                               "1, 0, 0, 0, 0, 0,"
                               "0, 1, 0, 0, 0, 0,"
                               "0, 0, 1, 0, 0, 0,"
                               "0, 0, 0,.5, 0, 0,"
                               "0, 0, 0, 0,.5, 0,"
                               "0, 0, 0, 0, 0,.5",
                               Patterns::List(Patterns::Double()),
                               "Viscosity-scaling tensor in Voigt notation.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Anisotropic<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);
      prm.enter_subsection ("Material model");
      {
        prm.enter_subsection ("Anisotropic");
        {
          const int size_voigt = (dim == 3 ? 6 : 3);
          const std::vector<double> tmp_tensor =
            Utilities::string_to_double(Utilities::split_string_list(prm.get ("Viscosity tensor")));
          Assert(tmp_tensor.size() == size_voigt*size_voigt,
                 ExcMessage("Constitutive voigt matrix must have 9 components in 2D, or 36 components in 3d"));

          std::vector<std::vector<double> > voigt_visc_tensor (size_voigt);
          for (unsigned int i=0; i<size_voigt; ++i)
            {
              voigt_visc_tensor[i].resize(size_voigt);
              for (unsigned int j=0; j<size_voigt; ++j)
                voigt_visc_tensor[i][j] = tmp_tensor[i*size_voigt+j];
            }

          // Voigt indices (For mapping back to real tensor)
          const unsigned int vi3d0[] = {0, 1, 2, 1, 0, 0};
          const unsigned int vi3d1[] = {0, 1, 2, 2, 2, 1};
          const unsigned int vi2d0[] = {0, 1, 0};
          const unsigned int vi2d1[] = {0, 1, 1};

          // Fill the constitutive tensor with values from the Voigt tensor
          for (unsigned int i=0; i<size_voigt; ++i)
            for (unsigned int j=0; j<size_voigt; ++j)
              if (dim == 2)
                C[vi2d0[i]][vi2d1[i]][vi2d0[j]][vi2d1[j]] = voigt_visc_tensor[i][j];
              else
                C[vi3d0[i]][vi3d1[i]][vi3d0[j]][vi3d1[j]] = voigt_visc_tensor[i][j];
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
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
  template class StokesCompressiblePreconditioner<dim>; \
  template class StokesIncompressibleTerms<dim>; \
  template class StokesCompressibleStrainRateViscosityTerm<dim>; \
  template class StokesReferenceDensityCompressibilityTerm<dim>; \
  template class StokesImplicitReferenceDensityCompressibilityTerm<dim>; \
  template class StokesIsentropicCompressionTerm<dim>; \
  template class StokesHydrostaticCompressionTerm<dim>; \
  template class StokesPressureRHSCompatibilityModification<dim>; \
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
    ASPECT_REGISTER_MATERIAL_MODEL(Anisotropic,
                                   "anisotropic",
                                   "A simple material model that is like the "
                                   "'Simple' model, but has a non-zero compressibility "
                                   "and always has a blob of anisotropic material in "
                                   "the center.")
  }
}
