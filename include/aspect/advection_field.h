/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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


#ifndef _aspect_advection_field_h
#define _aspect_advection_field_h

#include <aspect/global.h>
#include <aspect/parameters.h>

#include <deal.II/fe/fe_values_extractors.h>

namespace aspect
{
  using namespace dealii;

  /**
   * Forward declare the Introspection() class to avoid circular dependencies.
   */
  template <int dim> struct Introspection;

  /**
   * A structure that is used as an argument to functions that can work on
   * both the temperature and the compositional variables and that need to
   * be told which one of the two, as well as on which of the
   * compositional variables.
   */
  struct AdvectionField
  {
    /**
     * An enum indicating whether the identified variable is the
     * temperature or one of the compositional fields.
     */
    enum FieldType { temperature_field, compositional_field };

    /**
     * A variable indicating whether the identified variable is the
     * temperature or one of the compositional fields.
     */
    const FieldType    field_type;

    /**
     * A variable identifying which of the compositional fields is
     * selected. This variable is meaningless if the temperature is
     * selected.
     */
    const unsigned int compositional_variable;

    /**
     * Constructor.
     * @param field_type Determines whether this variable should select
     * the temperature field or a compositional field.
     * @param compositional_variable The number of the compositional field
     * if the first argument in fact chooses a compositional variable.
     * Meaningless if the first argument equals temperature.
     *
     * This function is implemented in
     * <code>source/simulator/helper_functions.cc</code>.
     */
    AdvectionField (const FieldType field_type,
                    const unsigned int compositional_variable = numbers::invalid_unsigned_int);

    /**
     * A static function that creates an object identifying the
     * temperature.
     *
     * This function is implemented in
     * <code>source/simulator/helper_functions.cc</code>.
     */
    static
    AdvectionField temperature ();

    /**
     * A static function that creates an object identifying given
     * compositional field.
     *
     * This function is implemented in
     * <code>source/simulator/helper_functions.cc</code>.
     */
    static
    AdvectionField composition (const unsigned int compositional_variable);

    /**
     * Return whether this object refers to the temperature field.
     */
    bool
    is_temperature () const;

    /**
     * Return whether this object refers to a field discretized by
     * discontinuous finite elements.
     */
    template<int dim>
    bool
    is_discontinuous (const Introspection<dim> &introspection) const;

    /**
     * Return the method that is used to solve the advection of this field
     * (i.e. 'fem_field', 'particles').
     */
    template<int dim>
    typename Parameters<dim>::AdvectionFieldMethod::Kind
    advection_method (const Introspection<dim> &introspection) const;

    /**
     * Look up the component index for this temperature or compositional
     * field. See Introspection::component_indices for more information.
     */
    template<int dim>
    unsigned int
    component_index(const Introspection<dim> &introspection) const;

    /**
     * Look up the block index for this temperature or compositional
     * field. See Introspection::block_indices for more information.
     */
    template<int dim>
    unsigned int
    block_index(const Introspection<dim> &introspection) const;

    /**
     * Look up the block index where the sparsity pattern for this field
     * is stored. This can be different than block_index() as several fields
     * can use the same pattern (typically in the first compositional field
     * if all fields are compatible). See Introspection::block_indices
     * for more information.
     */
    template<int dim>
    unsigned int
    sparsity_pattern_block_index(const Introspection<dim> &introspection) const;

    /**
     * Returns an index that runs from 0 (temperature field) to n (nth
     * compositional field), and uniquely identifies the current advection
     * field among the list of all advection fields. Can be used to index
     * vectors that contain entries for all advection fields.
     */
    unsigned int
    field_index() const;

    /**
     * Look up the base element within the larger composite finite element
     * we used for everything, for this temperature or compositional field
     * See Introspection::base_elements for more information.
     */
    template<int dim>
    unsigned int
    base_element(const Introspection<dim> &introspection) const;

    /**
     * Return the FEValues scalar extractor for this temperature
     * or compositional field.
     * This function is implemented in
     * <code>source/simulator/helper_functions.cc</code>.
     */
    template<int dim>
    FEValuesExtractors::Scalar
    scalar_extractor(const Introspection<dim> &introspection) const;

    /**
     * Look up the polynomial degree order for this temperature or compositional
     * field. See Introspection::polynomial_degree for more information.
     */
    template<int dim>
    unsigned int
    polynomial_degree(const Introspection<dim> &introspection) const;

    /**
     * Return a string that describes the field type and the compositional
     * variable number and name, if applicable.
     */
    template<int dim>
    std::string
    name(const Introspection<dim> &introspection) const;
  };
}


#endif
