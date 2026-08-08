/*
  Copyright (C) 2022 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/simple.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MaterialModel
  {
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
      if (in.current_cell.state() == IteratorState::valid && this->get_timestep_number() > 0)
        {
          const QGauss<dim> quadrature_formula (this->introspection().polynomial_degree.velocities +1);
          FEValues<dim> fe_values (this->get_mapping(),
                                   this->get_fe(),
                                   quadrature_formula,
                                   update_gradients);

          std::vector<Tensor<2,dim>> velocity_gradients (quadrature_formula.size(), Tensor<2,dim>());

          fe_values.reinit (in.current_cell);
          fe_values[this->introspection().extractors.velocities].get_function_gradients (this->get_solution(),
                                                                                         velocity_gradients);

          // Assign the strain components to the compositional fields reaction terms.
          // If there are too many fields, we simply fill only the first fields with the
          // existing strain rate tensor components.
          for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
            {
              // Convert the compositional fields into the tensor quantity they represent.
              const Tensor<2,dim> strain(make_array_view(&in.composition[q][0],
                                                         &in.composition[q][0] + Tensor<2,dim>::n_independent_components));

              // Compute the strain accumulated in this timestep.
              const Tensor<2,dim> strain_increment = this->get_timestep() * (velocity_gradients[q] * strain);

              // Output the strain increment component-wise to its respective compositional field's reaction terms.
              strain_increment.unroll(&out.reaction_terms[q][0],
                                      &out.reaction_terms[q][0] + Tensor<2,dim>::n_independent_components);
            }
        }
    }


    template <int dim>
    void
    FiniteStrain<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      Simple<dim>::parse_parameters (prm);

      AssertThrow(this->n_compositional_fields() >= (Tensor<2,dim>::n_independent_components),
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
