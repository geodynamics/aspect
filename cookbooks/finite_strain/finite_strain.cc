#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class FiniteStrain : public MaterialModel::Simple<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
        virtual void parse_parameters(ParameterHandler &prm);
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    FiniteStrain<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // First, we use the material descriptions of the 'simple' material model to fill all of the material
      // model outputs. Below, we will then overwrite selected properties (the reaction terms), which are
      // needed to track the finite strain.
      Simple<dim>::evaluate(in, out);

      // We need the velocity gradient for the finite strain (they are not included in material model inputs),
      // so we get them from the finite element.
      if (in.cell && this->get_timestep_number() > 0)
        {
          const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+1);
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   quadrature_formula,
                                   update_gradients);

          std::vector<Tensor<2,dim> > velocity_gradients (quadrature_formula.size(), Tensor<2,dim>());

          fe_values.reinit (*in.cell);
          fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_solution(),
                                                                                         velocity_gradients);

          // Assign the strain rate components to the compositional fields reaction terms.
          // We also have to rotate the accumulated strain from all the previous time steps with
          // the rotation part of the flow field of the current time step.
          // If there are too many fields, we simply fill only the first fields with the
          // existing strain rate tensor components.
          for (unsigned int q=0; q < in.position.size(); ++q)
            {
              // rotation tensor =
              //     asymmetric part of the displacement in this time step
              //                (= velocity gradient tensor * time step)
              //   + unit tensor
              const Tensor<2,dim> rotation = (velocity_gradients[q] - symmetrize(velocity_gradients[q])) * this->get_timestep()
                                             + unit_symmetric_tensor<dim>();

              SymmetricTensor<2,dim> accumulated_strain;
              for (unsigned int i=0; i<SymmetricTensor<2,dim>::n_independent_components; ++i)
                accumulated_strain[SymmetricTensor<2,dim>::unrolled_to_component_indices(i)] = in.composition[q][i];

              // the new strain is the rotated old strain plus the
              // strain of the current time step
              const SymmetricTensor<2,dim> rotated_strain = symmetrize(rotation * Tensor<2,dim>(accumulated_strain) * transpose(rotation)) + in.strain_rate[q] * this->get_timestep();

              for (unsigned int c=0; c<SymmetricTensor<2,dim>::n_independent_components; ++c)
                {
                  out.reaction_terms[q][c] = - in.composition[q][c]
                                             + rotated_strain[SymmetricTensor<2,dim>::unrolled_to_component_indices(c)];
                }
            }
        }
    }


    template <int dim>
    void
    FiniteStrain<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);

      AssertThrow(this->n_compositional_fields() >= (SymmetricTensor<2,dim>::n_independent_components),
                  ExcMessage("There must be at least as many compositional fields as independent components in the full "
                             "strain rate tensor."));
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(FiniteStrain,
                                   "finite strain",
                                   "A simple material model that is like the "
                                   "'Simple' model, but tracks the finite strain as compositional "
                                   "fields. The model assumes that the first 4 (in 2D) "
                                   " or 9 (in 3D) compositional fields contain the finite "
                                   "strain components. ")
  }
}
