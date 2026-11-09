// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/advection_field.h>
#include <aspect/introspection.h>
#include <aspect/parameters.h>


namespace aspect
{
  AdvectionField::AdvectionField (const FieldType field_type,
                                  const unsigned int compositional_variable)
    :
    field_type (field_type),
    compositional_variable (compositional_variable)
  {
    if (field_type == temperature_field)
      Assert (compositional_variable == numbers::invalid_unsigned_int,
              ExcMessage ("You can't specify a compositional variable if you "
                          "have in fact selected the temperature."));
  }



  AdvectionField
  AdvectionField::temperature ()
  {
    return AdvectionField(temperature_field);
  }



  AdvectionField
  AdvectionField::composition (const unsigned int compositional_variable)
  {
    return AdvectionField(compositional_field,
                          compositional_variable);
  }



  bool
  AdvectionField::is_temperature() const
  {
    return (field_type == temperature_field);
  }



  template <int dim>
  bool
  AdvectionField::is_discontinuous(const Introspection<dim> &introspection) const
  {
    if (field_type == temperature_field)
      return introspection.use_discontinuous_temperature_discretization;
    else if (field_type == compositional_field)
      return introspection.use_discontinuous_composition_discretization[compositional_variable];

    Assert (false, ExcInternalError());
    return false;
  }



  template <int dim>
  typename Parameters<dim>::AdvectionFieldMethod::Kind
  AdvectionField::advection_method(const Introspection<dim> &introspection) const
  {
    if (field_type == temperature_field)
      return introspection.temperature_method;
    else if (field_type == compositional_field)
      return introspection.compositional_field_methods[compositional_variable];

    Assert (false, ExcInternalError());
    return Parameters<dim>::AdvectionFieldMethod::fem_field;
  }



  template <int dim>
  unsigned int
  AdvectionField::block_index(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.block_indices.temperature;
    else
      return introspection.block_indices.compositional_fields[compositional_variable];
  }



  template <int dim>
  unsigned int
  AdvectionField::sparsity_pattern_block_index(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.block_indices.temperature;
    else
      return introspection.block_indices.compositional_field_sparsity_pattern[compositional_variable];
  }



  template <int dim>
  unsigned int
  AdvectionField::component_index(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.component_indices.temperature;
    else
      return introspection.component_indices.compositional_fields[compositional_variable];
  }



  unsigned int
  AdvectionField::field_index() const
  {
    if (this->is_temperature())
      return 0;
    else
      return compositional_variable + 1;
  }



  template <int dim>
  unsigned int
  AdvectionField::base_element(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.base_elements.temperature;
    else
      return introspection.base_elements.compositional_fields[compositional_variable];
  }



  template <int dim>
  FEValuesExtractors::Scalar
  AdvectionField::scalar_extractor(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.extractors.temperature;
    else
      {
        Assert(compositional_variable < introspection.n_compositional_fields,
               ExcMessage("Invalid AdvectionField."));
        return introspection.extractors.compositional_fields[compositional_variable];
      }
  }



  template <int dim>
  unsigned int
  AdvectionField::polynomial_degree(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return introspection.polynomial_degree.temperature;
    else
      return introspection.polynomial_degree.compositional_fields[compositional_variable];
  }



  template <int dim>
  std::string
  AdvectionField::name(const Introspection<dim> &introspection) const
  {
    if (this->is_temperature())
      return "temperature";
    else
      return "composition " + std::to_string(compositional_variable) + " (" + introspection.name_for_compositional_index(compositional_variable) + ")";
  }
}

// explicit instantiation of the functions we implement in this file
namespace aspect
{
#define INSTANTIATE(dim) \
  template bool AdvectionField::is_discontinuous (const Introspection<dim> &) const; \
  template typename Parameters<dim>::AdvectionFieldMethod::Kind AdvectionField::advection_method (const Introspection<dim> &) const; \
  template unsigned int AdvectionField::block_index (const Introspection<dim> &) const; \
  template unsigned int AdvectionField::sparsity_pattern_block_index (const Introspection<dim> &) const; \
  template unsigned int AdvectionField::component_index (const Introspection<dim> &) const; \
  template unsigned int AdvectionField::base_element(const Introspection<dim> &) const; \
  template FEValuesExtractors::Scalar AdvectionField::scalar_extractor (const Introspection<dim> &) const; \
  template unsigned int AdvectionField::polynomial_degree (const Introspection<dim> &) const; \
  template std::string AdvectionField::name (const Introspection<dim> &) const;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
