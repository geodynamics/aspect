/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>

namespace aspect
{
  namespace internal
  {
    /**
     * Return a filled ComponentIndices structure.
     */
    template <int dim>
    typename Introspection<dim>::ComponentIndices
    setup_component_indices (FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::ComponentIndices ci;
      for (unsigned int i=0; i<dim; ++i)
        ci.velocities[i] = fevs.variable("velocity").first_component_index+i;

      ci.pressure = fevs.variable("pressure").first_component_index;
      ci.temperature = fevs.variable("temperature").first_component_index;

      {
        const auto composition_variables = fevs.variables_with_name("compositions");
        for (const auto &fe : composition_variables)
          {
            for (unsigned int i=0; i<fe->multiplicity; ++i)
              {
                ci.compositional_fields.push_back(fe->first_component_index+i);
              }
          }
      }

      return ci;
    }

    /**
     * return filled BlockIndices structure.
     */
    template <int dim>
    typename Introspection<dim>::BlockIndices
    setup_blocks (FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::BlockIndices b;
      b.velocities = fevs.variable("velocity").block_index;
      b.pressure = fevs.variable("pressure").block_index;
      b.temperature = fevs.variable("temperature").block_index;

      {
        const auto variables = fevs.variables_with_name("compositions");
        for (const auto &fe : variables)
          {
            for (unsigned int i=0; i<fe->multiplicity; ++i)
              {
                b.compositional_fields.push_back(fe->block_index + i);
                b.compositional_field_sparsity_pattern.push_back(fe->block_index);
              }
          }
      }
      return b;
    }



    template <int dim>
    typename Introspection<dim>::BaseElements
    setup_base_elements (FEVariableCollection<dim> &fevs)
    {
      typename Introspection<dim>::BaseElements base_elements;

      base_elements.velocities = fevs.variable("velocity").base_index;
      base_elements.pressure = fevs.variable("pressure").base_index;
      base_elements.temperature = fevs.variable("temperature").base_index;

      {
        // TODO: We would ideally reuse base elements also when they are
        // not consecutively encountered. For now, just ignore that fact.
        const auto variables = fevs.variables_with_name("compositions");
        for (const auto &fe : variables)
          {
            for (unsigned int i=0; i<fe->multiplicity; ++i)
              base_elements.compositional_fields.push_back(fe->base_index);
          }
      }




      return base_elements;
    }



    template <int dim>
    typename Introspection<dim>::PolynomialDegree
    setup_polynomial_degree (const Parameters<dim> &parameters)
    {
      typename Introspection<dim>::PolynomialDegree polynomial_degree;

      polynomial_degree.velocities = parameters.stokes_velocity_degree;
      polynomial_degree.temperature = parameters.temperature_degree;
      polynomial_degree.compositional_fields = parameters.composition_degrees;
      polynomial_degree.max_compositional_field = parameters.max_composition_degree;
      polynomial_degree.max_degree = std::max({parameters.stokes_velocity_degree,
                                               parameters.temperature_degree,
                                               (parameters.n_compositional_fields>0 ? parameters.max_composition_degree : 0u)
                                              }
                                             );

      return polynomial_degree;
    }



    template <int dim>
    typename Introspection<dim>::Quadratures
    setup_quadratures (const Parameters<dim> &parameters,
                       const ReferenceCell reference_cell)
    {
      typename Introspection<dim>::Quadratures quadratures;

      quadratures.velocities = reference_cell.get_gauss_type_quadrature<dim>(parameters.stokes_velocity_degree+1);
      quadratures.pressure = reference_cell.get_gauss_type_quadrature<dim>(parameters.use_equal_order_interpolation_for_stokes
                                                                           ?
                                                                           parameters.stokes_velocity_degree+1
                                                                           :
                                                                           parameters.stokes_velocity_degree);
      quadratures.temperature = reference_cell.get_gauss_type_quadrature<dim>(parameters.temperature_degree+1);
      quadratures.compositional_field_max = reference_cell.get_gauss_type_quadrature<dim>(parameters.max_composition_degree+1);
      for (const auto &degree: parameters.composition_degrees)
        quadratures.compositional_fields.emplace_back(reference_cell.get_gauss_type_quadrature<dim>(degree+1));
      quadratures.system = reference_cell.get_gauss_type_quadrature<dim>(std::max({parameters.stokes_velocity_degree,
                                                                                   parameters.temperature_degree,
                                                                                   parameters.max_composition_degree
                                                                                  }) + 1);

      return quadratures;
    }



    template <int dim>
    typename Introspection<dim>::FaceQuadratures
    setup_face_quadratures (const Parameters<dim> &parameters,
                            const ReferenceCell reference_cell)
    {
      typename Introspection<dim>::FaceQuadratures quadratures;

      quadratures.velocities = reference_cell.face_reference_cell(0).get_gauss_type_quadrature<dim-1>(parameters.stokes_velocity_degree+1);
      quadratures.pressure = reference_cell.face_reference_cell(0).get_gauss_type_quadrature<dim-1>(parameters.use_equal_order_interpolation_for_stokes
                             ?
                             parameters.stokes_velocity_degree+1
                             :
                             parameters.stokes_velocity_degree);
      quadratures.temperature = reference_cell.face_reference_cell(0).get_gauss_type_quadrature<dim-1>(parameters.temperature_degree+1);
      quadratures.compositional_fields = reference_cell.face_reference_cell(0).get_gauss_type_quadrature<dim-1>(parameters.max_composition_degree+1);
      quadratures.system = reference_cell.face_reference_cell(0).get_gauss_type_quadrature<dim-1>(
                             std::max({parameters.stokes_velocity_degree, parameters.temperature_degree, parameters.max_composition_degree}) + 1);

      return quadratures;
    }



    template <int dim>
    std::shared_ptr<FiniteElement<dim>>
    new_FE_Q_or_DGP(const bool discontinuous,
                    const unsigned int degree)
    {
      if (discontinuous)
        return std::make_shared<FE_DGP<dim>>(degree);
      else
        return std::make_shared<FE_Q<dim>>(degree);
    }



    template <int dim>
    std::shared_ptr<FiniteElement<dim>>
    new_FE_Q_or_DGQ(const bool discontinuous,
                    const unsigned int degree)
    {
      if (discontinuous)
        return std::make_shared<FE_DGQ<dim>>(degree);
      else
        return std::make_shared<FE_Q<dim>>(degree);
    }
  }




  template <int dim>
  std::vector<VariableDeclaration<dim>>
  construct_default_variables (const Parameters<dim> &parameters)
  {
    std::vector<VariableDeclaration<dim>> variables;

    const unsigned int n_velocity_blocks = parameters.use_direct_stokes_solver ? 0 : 1;
    variables.push_back(
      VariableDeclaration<dim>("velocity",
                               std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree),
                               dim,
                               n_velocity_blocks));

    if (parameters.use_equal_order_interpolation_for_stokes == false)
      variables.push_back(
        VariableDeclaration<dim>(
          "pressure",
          internal::new_FE_Q_or_DGP<dim>(parameters.use_locally_conservative_discretization,
                                         parameters.stokes_velocity_degree-1),
          1,
          1));
    else
      variables.push_back(
        VariableDeclaration<dim>(
          "pressure",
          std::make_shared<FE_Q<dim>>(parameters.stokes_velocity_degree),
          1,
          1));


    variables.push_back(
      VariableDeclaration<dim>(
        "temperature",
        internal::new_FE_Q_or_DGQ<dim>(parameters.use_discontinuous_temperature_discretization,
                                       parameters.temperature_degree),
        1,
        1));

    // We can not create a single composition with multiplicity n (this assumes all have the same FiniteElement) and
    // we also don't want to create n different compositional fields for performance reasons. Instead, we combine
    // consecutive compositions of the same type into the same base element. For example, if the user requests
    // "Q2,Q2,DGQ1,Q2", we actually create "Q2^2, DGQ1, Q2":
    for (unsigned int c=0; c<parameters.n_compositional_fields; ++c)
      {
        if (c>0
            && parameters.use_discontinuous_composition_discretization[c] == parameters.use_discontinuous_composition_discretization[c-1]
            && parameters.composition_degrees[c] == parameters.composition_degrees[c-1])
          {
            // reuse last one because it is the same
            variables.back().multiplicity += 1;
            variables.back().n_blocks += 1;
          }
        else
          {
            std::shared_ptr<FiniteElement<dim>> fe = internal::new_FE_Q_or_DGQ<dim>(parameters.use_discontinuous_composition_discretization[c],
                                                                                     parameters.composition_degrees[c]);
            variables.push_back(VariableDeclaration<dim>("compositions", fe, 1, 1));
          }
      }

    return variables;
  }


  template <int dim>
  Introspection<dim>::Introspection(const std::vector<VariableDeclaration<dim>> &variable_definition,
                                    const Parameters<dim> &parameters
                                   )
    :
    FEVariableCollection<dim>(variable_definition),
    n_components (FEVariableCollection<dim>::n_components()),
    n_compositional_fields (parameters.n_compositional_fields),
    use_discontinuous_temperature_discretization (parameters.use_discontinuous_temperature_discretization),
    use_discontinuous_composition_discretization (parameters.use_discontinuous_composition_discretization),
    component_indices (internal::setup_component_indices<dim>(*this)),
    n_blocks (FEVariableCollection<dim>::n_blocks()),
    block_indices (internal::setup_blocks<dim>(*this)),
    extractors (component_indices),
    base_elements (internal::setup_base_elements<dim>(*this)),
    polynomial_degree (internal::setup_polynomial_degree<dim>(parameters)),
    quadratures (internal::setup_quadratures<dim>(parameters, ReferenceCells::get_hypercube<dim>())),
    face_quadratures (internal::setup_face_quadratures<dim>(parameters, ReferenceCells::get_hypercube<dim>())),
    component_masks (*this, component_indices),
    system_dofs_per_block (n_blocks),
    temperature_method(parameters.temperature_method),
    compositional_field_methods(parameters.compositional_field_methods),
    composition_names(parameters.names_of_compositional_fields),
    composition_descriptions(parameters.composition_descriptions),
    composition_names_for_type(CompositionalFieldDescription::n_types),
    composition_indices_for_type(CompositionalFieldDescription::n_types)
  {
    // Set up the indices for the different types of compositional fields
    for (unsigned int c=0; c<composition_descriptions.size(); ++c)
      {
        composition_indices_for_type[composition_descriptions[c].type].push_back(c);
        composition_names_for_type[composition_descriptions[c].type].push_back(composition_names[c]);
      }

    // Fill composition_base_element_indices
    {
      if (this->n_compositional_fields > 0)
        {
          // We are assigning base elements in order, so the first compositional field
          // gives us the first base element index. Then we find the largest index
          // in the vector. This is necessary, because the fields could have type A,B,A.
          const unsigned int first = this->base_elements.compositional_fields[0];
          const unsigned int last = *std::max_element(this->base_elements.compositional_fields.begin(),
                                                      this->base_elements.compositional_fields.end());

          composition_base_element_indices.resize(last-first+1);
          std::iota(composition_base_element_indices.begin(), composition_base_element_indices.end(), first);
        }
    }

// Fill compositional_field_indices_with_base_element
    {
      for (const auto base_element_index : composition_base_element_indices)
        {
          std::vector<unsigned int> result;

          unsigned int idx = 0;
          for (const auto base_idx : this->base_elements.compositional_fields)
            {
              if (base_idx == base_element_index)
                result.emplace_back(idx);
              ++idx;
            }

          Assert(result.size() > 0, ExcInternalError("There should be at least one compositional field for a valid base element."));
          compositional_field_indices_with_base_element[base_element_index] = result;
        }
    }

  }



  template <int dim>
  Introspection<dim>::~Introspection ()
    = default;



  namespace
  {
    std::vector<FEValuesExtractors::Scalar>
    make_extractor_sequence (const std::vector<unsigned int> &compositional_fields)
    {
      std::vector<FEValuesExtractors::Scalar> x;
      x.reserve(compositional_fields.size());
      for (const unsigned int compositional_field : compositional_fields)
        x.emplace_back(compositional_field);
      return x;
    }
  }



  template <int dim>
  Introspection<dim>::Extractors::Extractors (const Introspection<dim>::ComponentIndices &component_indices)
    :
    velocities (component_indices.velocities[0]),
    pressure (component_indices.pressure),
    temperature (component_indices.temperature),
    compositional_fields (make_extractor_sequence (component_indices.compositional_fields))
  {}



  namespace
  {
    template <int dim>
    std::vector<ComponentMask>
    make_composition_component_mask_sequence (const FEVariableCollection<dim> &fevs, const typename Introspection<dim>::ComponentIndices &indices)
    {
      std::vector<ComponentMask> result;
      for (const unsigned int idx : indices.compositional_fields)
        {
          result.push_back(ComponentMask(fevs.n_components(), false));
          result.back().set(idx, true);
        }
      return result;
    }

    template <int dim>
    ComponentMask
    make_composition_component_mask (const FEVariableCollection<dim> &fevs, const std::vector<ComponentMask> &compositional_fields)
    {
      ComponentMask result(fevs.n_components(), false);
      for (const auto &mask : compositional_fields)
        result = result | mask;
      return result;
    }
  }



  template <int dim>
  Introspection<dim>::ComponentMasks::ComponentMasks (const FEVariableCollection<dim> &fevs, const Introspection<dim>::ComponentIndices &indices)
    :
    velocities (fevs.variable("velocity").component_mask),
    pressure (fevs.variable("pressure").component_mask),
    temperature (fevs.variable("temperature").component_mask),
    compositional_fields (make_composition_component_mask_sequence (fevs, indices)),
    compositions(make_composition_component_mask(fevs, compositional_fields))
  {}



  template <int dim>
  const std::vector<unsigned int> &
  Introspection<dim>::get_composition_base_element_indices() const
  {
    return composition_base_element_indices;
  }



  template <int dim>
  const std::vector<unsigned int> &
  Introspection<dim>::get_compositional_field_indices_with_base_element(const unsigned int base_element_index) const
  {
    Assert(compositional_field_indices_with_base_element.find(base_element_index)
           != compositional_field_indices_with_base_element.end(),
           ExcMessage("Invalid base_element_index specified."));
    return compositional_field_indices_with_base_element.find(base_element_index)->second;
  }



  template <int dim>
  unsigned int
  Introspection<dim>::compositional_index_for_name (const std::string &name) const
  {
    const std::vector<std::string>::const_iterator
    it = std::find(composition_names.begin(), composition_names.end(), name);
    AssertThrow (it != composition_names.end(),
                 ExcMessage ("The compositional field " + name +
                             " you asked for is not used in the simulation."));

    return (it - composition_names.begin());
  }



  template <int dim>
  std::string
  Introspection<dim>::name_for_compositional_index (const unsigned int index) const
  {
    // make sure that what we get here is really an index of one of the compositional fields
    AssertIndexRange(index,composition_names.size());
    return composition_names[index];
  }



  template <int dim>
  const std::vector<std::string> &
  Introspection<dim>::get_composition_names () const
  {
    // Simply return the full list of composition names
    return composition_names;
  }



  template <int dim>
  const std::vector<CompositionalFieldDescription> &
  Introspection<dim>::get_composition_descriptions () const
  {
    return composition_descriptions;
  }



  template <int dim>
  const std::vector<std::string> &
  Introspection<dim>::chemical_composition_field_names () const
  {
    return get_names_for_fields_of_type(CompositionalFieldDescription::chemical_composition);
  }



  template <int dim>
  const std::vector<unsigned int> &
  Introspection<dim>::chemical_composition_field_indices () const
  {
    return get_indices_for_fields_of_type(CompositionalFieldDescription::chemical_composition);
  }



  template <int dim>
  unsigned int
  Introspection<dim>::n_chemical_composition_fields () const
  {
    return get_number_of_fields_of_type(CompositionalFieldDescription::Type::chemical_composition);
  }



  template <int dim>
  bool
  Introspection<dim>::composition_type_exists (const CompositionalFieldDescription::Type &type) const
  {
    Assert(type < composition_indices_for_type.size(), ExcInternalError());

    return composition_indices_for_type[type].size() > 0;
  }



  template <int dim>
  unsigned int
  Introspection<dim>::find_composition_type (const typename CompositionalFieldDescription::Type &type) const
  {
    Assert(type < composition_indices_for_type.size(), ExcInternalError());

    if (composition_indices_for_type[type].size() > 0)
      return composition_indices_for_type[type][0];

    return composition_descriptions.size();
  }



  template <int dim>
  bool
  Introspection<dim>::compositional_name_exists (const std::string &name) const
  {
    return (std::find(composition_names.begin(), composition_names.end(), name) != composition_names.end()
            ?
            true
            :
            false);
  }



  template <int dim>
  const std::vector<unsigned int> &
  Introspection<dim>::get_indices_for_fields_of_type (const CompositionalFieldDescription::Type &type) const
  {
    Assert(type < composition_indices_for_type.size(), ExcInternalError());
    return composition_indices_for_type[type];
  }



  template <int dim>
  const std::vector<std::string> &
  Introspection<dim>::get_names_for_fields_of_type (const CompositionalFieldDescription::Type &type) const
  {
    Assert(type < composition_names_for_type.size(), ExcInternalError());
    return composition_names_for_type[type];
  }



  template <int dim>
  unsigned int
  Introspection<dim>::get_number_of_fields_of_type (const CompositionalFieldDescription::Type &type) const
  {
    Assert(type < composition_indices_for_type.size(), ExcInternalError());
    return composition_indices_for_type[type].size();
  }



  template <int dim>
  bool
  Introspection<dim>::is_stokes_component (const unsigned int component_index) const
  {
    if (component_index == component_indices.pressure)
      return true;

    for (unsigned int i=0; i<dim; ++i)
      if (component_index == component_indices.velocities[i])
        return true;

    return false;
  }



  template <int dim>
  bool
  Introspection<dim>::is_composition_component (const unsigned int component_index) const
  {
    // All compositions live at the end. Just to be sure, there are no other components
    // in our system after compositional fields, right?
    Assert(component_indices.compositional_fields[0] > component_indices.temperature
           && component_indices.compositional_fields.back() == n_components-1, ExcInternalError());

    if (component_index >= component_indices.compositional_fields[0])
      return true;
    else
      return false;
  }

}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Introspection<dim>; \
  template \
  std::vector<VariableDeclaration<dim>> \
  construct_default_variables (const Parameters<dim> &parameters);

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
