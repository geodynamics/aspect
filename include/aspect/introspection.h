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


#ifndef _aspect_introspection_h
#define _aspect_introspection_h

#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe.h>

#include <aspect/fe_variable_collection.h>
#include <aspect/parameters.h>

#include <map>

namespace aspect
{
  /**
   * Helper function to construct the default list of variables to use
   * based on the given set of @p parameters.
   */
  template <int dim>
  std::vector<VariableDeclaration<dim>>
  construct_default_variables (const Parameters<dim> &parameters);


  /**
   * A data structure containing a description of each compositional field.
   * At present, this structure only includes the field type
   * (i.e., whether it is of type chemical composition, porosity, etc.).
   */
  struct CompositionalFieldDescription
  {
    /**
     * This enum lists available compositional field types.
     */
    enum Type
    {
      chemical_composition = 0,
      stress = 1,
      strain = 2,
      grain_size = 3,
      porosity = 4,
      density = 5,
      entropy = 6,
      generic = 7,
      unspecified = 8
    } type;

    /**
     * The number of different types defined in Type.
     */
    constexpr static unsigned int n_types = 9;

    /**
     * This function translates an input string into the
     * available enum options for the type of compositional field.
     */
    static
    Type
    parse_type(const std::string &input)
    {
      if (input == "chemical composition")
        return CompositionalFieldDescription::chemical_composition;
      else if (input == "stress")
        return CompositionalFieldDescription::stress;
      else if (input == "strain")
        return CompositionalFieldDescription::strain;
      else if (input == "grain size")
        return CompositionalFieldDescription::grain_size;
      else if (input == "porosity")
        return CompositionalFieldDescription::porosity;
      else if (input == "density")
        return CompositionalFieldDescription::density;
      else if (input == "entropy")
        return CompositionalFieldDescription::entropy;
      else if (input == "generic")
        return CompositionalFieldDescription::generic;
      else if (input == "unspecified")
        return CompositionalFieldDescription::unspecified;
      else
        AssertThrow(false, ExcNotImplemented());

      return CompositionalFieldDescription::Type();
    }
  };

  /**
   * The introspection class provides information about the simulation as a
   * whole. In particular, it provides symbolic names for things like the
   * velocity, pressure and other variables, along with their corresponding
   * vector and scalar FEValues extractors, component masks, etc.
   *
   * The purpose of this class is primarily to provide these symbolic names so
   * that we do not have to use implicit knowledge about the ordering of
   * variables (e.g., earlier versions of ASPECT had many places where we
   * built a scalar FEValues extractor at component 'dim' since that is where
   * we knew that the pressure lies in the finite element; this kind of
   * implicit knowledge is no longer necessary with the Introspection class).
   * The Introspection class is used both by the Simulator class itself, but
   * is also exported to plugins via the SimulatorAccess class.
   *
   * @ingroup Simulator
   */
  template <int dim>
  struct Introspection: public FEVariableCollection<dim>
  {
    public:
      /**
       * Constructor.
       */
      Introspection (const std::vector<VariableDeclaration<dim>> &variables,
                     const Parameters<dim> &parameters);

      /**
       * Destructor.
       */
      ~Introspection ();



      /**
       * @name Things that are independent of the current mesh
       * @{
       */

      /**
       * The number of vector components used by the finite element
       * description of this problem. It equals $d+2+n_c$ where $d$ is the
       * dimension and equals the number of velocity components, and $n_c$ is
       * the number of advected (compositional) fields. The remaining
       * components are the scalar pressure and temperature fields.
       */
      const unsigned int n_components;

      /**
       * The number of compositional fields.
       */
      const unsigned int n_compositional_fields;

      /**
       * A variable that holds whether the temperature field should use a
       * discontinuous discretization.
       */
      const bool use_discontinuous_temperature_discretization;

      /**
       * A variable that holds whether the composition field(s) should use a
       * discontinuous discretization.
       */
      const std::vector<bool> use_discontinuous_composition_discretization;

      /**
       * A structure that enumerates the vector components of the finite
       * element that correspond to each of the variables in this problem.
       */
      struct ComponentIndices
      {
        std::array<unsigned int, dim> velocities;
        unsigned int                  pressure;
        unsigned int                  temperature;
        std::vector<unsigned int>     compositional_fields;
      };
      /**
       * A variable that enumerates the vector components of the finite
       * element that correspond to each of the variables in this problem.
       */
      const ComponentIndices component_indices;

      /**
       * The number of vector blocks. This equals $3+n_c$ where, in comparison
       * to the n_components field, the velocity components form a single
       * block.
       */
      const unsigned int n_blocks;

      /**
       * A structure that enumerates the vector blocks of the finite element
       * that correspond to each of the variables in this problem.
       */
      struct BlockIndices
      {
        unsigned int       velocities;
        unsigned int       pressure;
        unsigned int       temperature;
        std::vector<unsigned int> compositional_fields;

        /**
         * This variable contains the block for each compositional field
         * where the matrix/sparsity pattern is copied from when we need to
         * (temporarily) create a matrix. This way, we only need to store a
         * single sparsity pattern and reuse it for all compositional fields
         * (assuming they have an identical FiniteElement).
         */
        std::vector<unsigned int> compositional_field_sparsity_pattern;
      };

      /**
       * A variable that enumerates the vector blocks of the finite element
       * that correspond to each of the variables in this problem.
       */
      const BlockIndices block_indices;

      /**
       * A structure that contains FEValues extractors for every block of the
       * finite element used in the overall description.
       */
      struct Extractors
      {
        Extractors (const ComponentIndices &component_indices);

        const FEValuesExtractors::Vector              velocities;
        const FEValuesExtractors::Scalar              pressure;
        const FEValuesExtractors::Scalar              temperature;
        const std::vector<FEValuesExtractors::Scalar> compositional_fields;
      };

      /**
       * A variable that contains extractors for every block of the finite
       * element used in the overall description.
       */
      const Extractors extractors;

      /**
       * A structure that enumerates the base elements of the finite element
       * that correspond to each of the variables in this problem.
       *
       * The indices here can be used to access the dealii::FiniteElement
       * that describes the given variable. We support different finite
       * elements for compositional fields, but we try to reuse the same
       * element if possible.
       */
      struct BaseElements
      {
        unsigned int              velocities;
        unsigned int              pressure;
        unsigned int              temperature;
        std::vector<unsigned int> compositional_fields;
      };

      /**
       * A variable that enumerates the base elements of the finite element
       * that correspond to each of the variables in this problem.
       */
      const BaseElements base_elements;

      /**
       * A structure that contains the polynomial degree of the finite element
       * that correspond to each of the variables in this problem.
       *
       * If there are compositional fields, they are all discretized with the
       * same polynomial degree and, consequently, we only need a single integer.
       */
      struct PolynomialDegree
      {
        unsigned int              max_degree;
        unsigned int              velocities;
        unsigned int              temperature;
        std::vector<unsigned int> compositional_fields;
        unsigned int              max_compositional_field;
      };

      /**
       * A variable that enumerates the polynomial degree of the finite element
       * that correspond to each of the variables in this problem.
       */
      const PolynomialDegree polynomial_degree;

      /**
       * A structure that contains appropriate quadrature formulas for the
       * finite elements that correspond to each of the variables in this problem,
       * as well as for the complete system (the `system` variable).
       *
       * If there are compositional fields, they are all discretized with the
       * same polynomial degree and, consequently, we only need a single formula.
       *
       * The quadrature formulas provided here are chosen such that the
       * compute integrals with sufficient accuracy. For hypercube cells,
       * this means in particular that they are of Gauss type with a
       * number of Gauss points per coordinate direction that is one larger than
       * the polynomial degree of the finite element per direction. For
       * example, when using quadratic elements for the velocity, the
       * corresponding quadrature formula will have three Gauss points per
       * direction. If the mesh is based on triangles or tetrahedra, the
       * quadrature formula is not of tensor-product Gauss type, but the
       * corresponding analog for simplex cells.
       */
      struct Quadratures
      {
        Quadrature<dim>       velocities;
        Quadrature<dim>       pressure;
        Quadrature<dim>       temperature;
        Quadrature<dim>       compositional_field_max;
        std::vector<Quadrature<dim>> compositional_fields;
        Quadrature<dim>       system;
      };

      /**
       * A variable that enumerates the polynomial degree of the finite element
       * that correspond to each of the variables in this problem.
       */
      const Quadratures quadratures;

      /**
       * A structure that contains appropriate face quadrature formulas for the
       * finite elements that correspond to each of the variables in this problem,
       * as well as for the complete system (the `system` variable).
       *
       * This structure corresponds to the Quadratures structure above, but for
       * face integration.
       */
      struct FaceQuadratures
      {
        Quadrature<dim-1>       velocities;
        Quadrature<dim-1>       pressure;
        Quadrature<dim-1>       temperature;
        Quadrature<dim-1>       compositional_fields;
        Quadrature<dim-1>       system;
      };

      /**
       * A variable that enumerates the polynomial degree of the finite element
       * that correspond to each of the variables in this problem.
       */
      const FaceQuadratures face_quadratures;

      /**
       * A structure that contains component masks for each of the variables
       * in this problem. Component masks are a deal.II concept, see the
       * deal.II glossary.
       */
      struct ComponentMasks
      {
        ComponentMasks (const FEVariableCollection<dim> &fevs, const Introspection<dim>::ComponentIndices &indices);

        /**
         * The component mask for all velocity components.
         */
        ComponentMask              velocities;

        /**
         * The component mask for the pressure component.
         */
        ComponentMask              pressure;

        /**
         * The component mask for the temperature component.
         */
        ComponentMask              temperature;

        /**
         * The component mask for each individual compositional field.
         * The size of this vector is equal to the number of compositional fields.
         * Each entry is a component mask that selects the component
         * that corresponds to the respective compositional field.
         */
        std::vector<ComponentMask> compositional_fields;

        /**
         * The component mask for all composition components.
         * This mask selects all compositional fields.
         */
        ComponentMask              compositions;
      };

      /**
       * A variable that contains component masks for each of the variables in
       * this problem. Component masks are a deal.II concept, see the deal.II
       * glossary.
       */
      const ComponentMasks component_masks;

      /**
       * @}
       */

      /**
       * @name Things that depend on the current mesh
       * @{
       */
      /**
       * A variable that describes how many of the degrees of freedom on the
       * current mesh belong to each of the n_blocks blocks of the finite
       * element.
       */
      std::vector<types::global_dof_index> system_dofs_per_block;

      /**
       * A structure that contains index sets describing which of the globally
       * enumerated degrees of freedom are owned by or are relevant to the
       * current processor in a parallel computation.
       */
      struct IndexSets
      {
        /**
         * An index set that indicates which (among all) degrees of freedom
         * are relevant to the current processor. See the deal.II
         * documentation for the definition of the term "locally relevant
         * degrees of freedom".
         */
        IndexSet system_relevant_set;

        /**
         * A collection of index sets that for each of the vector blocks of
         * this finite element represents the global indices of the degrees of
         * freedom owned by this processor. The n_blocks elements of this
         * array form a mutually exclusive decomposition of the index set
         * containing all locally owned degrees of freedom.
         */
        std::vector<IndexSet> system_partitioning;

        /**
         * A collection of index sets that for each of the vector blocks of
         * this finite element represents the global indices of the degrees of
         * freedom are relevant to this processor. The n_blocks elements of
         * this array form a mutually exclusive decomposition of the index set
         * containing all locally relevant degrees of freedom, i.e., of the
         * system_relevant_set index set.
         */
        std::vector<IndexSet> system_relevant_partitioning;

        /**
         * A collection of index sets for each vector block of the Stokes
         * system (velocity and pressure). This variable contains the first
         * two elements of system_partitioning.
         */
        std::vector<IndexSet> stokes_partitioning;

        /**
         * Pressure unknowns that are locally owned. This IndexSet is needed
         * if velocity and pressure end up in the same block and is used for
         * pressure scaling and in make_pressure_rhs_compatible(). If melt
         * transport is enabled, this field is unused and not filled.
         */
        IndexSet locally_owned_pressure_dofs;

        /**
         * Fluid and compaction pressure unknowns that are locally owned. Only
         * valid if melt transport is enabled.
         */
        IndexSet locally_owned_melt_pressure_dofs;

        /**
         * Fluid pressure unknowns that are locally owned. Only valid if melt
         * transport is enabled.
         */
        IndexSet locally_owned_fluid_pressure_dofs;
      };

      /**
       * A variable that contains index sets describing which of the globally
       * enumerated degrees of freedom are owned by the current processor in a
       * parallel computation.
       */
      IndexSets index_sets;

      /**
       * A variable that contains the field method for the temperature field
       * and is used to determine how to solve it when solving a timestep.
       */
      typename Parameters<dim>::AdvectionFieldMethod::Kind temperature_method;

      /**
       * A vector that contains a field method for every compositional
       * field and is used to determine how to solve a particular field when
       * solving a timestep.
       */
      std::vector<typename Parameters<dim>::AdvectionFieldMethod::Kind> compositional_field_methods;

      /**
       * @}
       */

      /**
       * Return a vector that contains the base element indices of the deal.II FiniteElement
       * that are used for compositional fields. Note that compositional fields can share the
       * same base element, so this vector can (and usually will) be smaller than the number
       * of compositional fields. The function get_compositional_field_indices_with_base_element()
       * can be used to translate from base element index to all compositional field indices
       * using the specified base element.
       * If several fields are the finite element type (same degree and both continuous or both
       * discontinuous), they share base elements. If you have no compositional fields, the
       * vector returned has length 0. If all compositional fields have the same finite element
       * space, the length is 1.
       */
      const std::vector<unsigned int> &
      get_composition_base_element_indices() const;

      /**
       * Return a vector with all compositional field indices that belong to a given
       * base element index as returned by get_composition_base_element_indices().
       * The indices returned are therefore between 0 and n_compositional_fields-1.
       * If you have a single compositional field, this function returns {0} when passing
       * in base_element_index=0.
       */
      const std::vector<unsigned int> &
      get_compositional_field_indices_with_base_element(const unsigned int base_element_index) const;

      /**
       * A function that gets the name of a compositional field as an input
       * parameter and returns its index. If the name is not found, an
       * exception is thrown.
       *
       * @param name The name of compositional field (as specified in the
       * input file)
       */
      unsigned int
      compositional_index_for_name (const std::string &name) const;

      /**
       * A function that gets the index of a compositional field as an input
       * parameter and returns its name.
       *
       * @param index The index of compositional field
       */
      std::string
      name_for_compositional_index (const unsigned int index) const;

      /**
       * A function that returns the full list of compositional field names.
       */
      const std::vector<std::string> &
      get_composition_names () const;

      /**
       * A function that returns the full vector of compositional
       * field descriptions.
       */
      const std::vector<CompositionalFieldDescription> &
      get_composition_descriptions () const;

      /**
       * A function that returns the names of
       * compositional fields that correspond to chemical compositions.
       *
       * This function is shorthand for
       * get_names_for_fields_of_type(CompositionalFieldDescription::chemical_composition).
       */
      const std::vector<std::string> &
      chemical_composition_field_names () const;

      /**
       * A function that returns the indices of
       * compositional fields that correspond to chemical compositions.
       *
       * This function is shorthand for
       * get_indices_for_fields_of_type(CompositionalFieldDescription::chemical_composition).
       */
      const std::vector<unsigned int> &
      chemical_composition_field_indices () const;

      /**
       * A function that returns the number of
       * compositional fields that correspond to chemical compositions.
       *
       * This function is shorthand for
       * get_number_of_fields_of_type(CompositionalFieldDescription::chemical_composition).
       */
      unsigned int
      n_chemical_composition_fields () const;

      /**
       * A function that gets the type of a compositional field as an input
       * parameter and returns if any compositional field of that type is
       * used in this simulation.
       *
       * @param type The type of compositional field (as specified in the
       * input file)
       */
      bool
      composition_type_exists (const CompositionalFieldDescription::Type &type) const;

      /**
       * A function that gets the type of a compositional field as an input
       * parameter and returns the index of the first compositional field of
       * this type used in this simulation. If no such field is found, the
       * function returns the number of compositional fields.
       *
       * @param type The type of compositional field (as specified in the
       * input file)
       */
      unsigned int
      find_composition_type (const CompositionalFieldDescription::Type &type) const;

      /**
       * A function that gets the name of a compositional field as an input
       * parameter and returns if the compositional field is used in this
       * simulation.
       *
       * @param name The name of compositional field (as specified in the
       * input file)
       */
      bool
      compositional_name_exists (const std::string &name) const;

      /**
       * Get the indices of the compositional fields which are of a
       * particular type (chemical composition, porosity, etc.).
       */
      const std::vector<unsigned int> &
      get_indices_for_fields_of_type (const CompositionalFieldDescription::Type &type) const;

      /**
       * Get the names of the compositional fields which are of a
       * particular type (chemical composition, porosity, etc.).
       */
      const std::vector<std::string> &
      get_names_for_fields_of_type (const CompositionalFieldDescription::Type &type) const;

      /**
       * Get the number of compositional fields which are of a
       * particular type (chemical composition, porosity, etc.).
       */
      unsigned int
      get_number_of_fields_of_type (const CompositionalFieldDescription::Type &type) const;

      /**
       * A function that gets a component index as an input
       * parameter and returns if the component is one of the stokes system
       * (i.e. if it is the pressure or one of the velocity components).
       *
       * @param component_index The component index to check.
       */
      bool
      is_stokes_component (const unsigned int component_index) const;

      /**
       * A function that gets a component index as an input
       * parameter and returns if the component is one of the
       * compositional fields.
       *
       * @param component_index The component index to check.
       */
      bool
      is_composition_component (const unsigned int component_index) const;

    private:
      /**
       * A vector that stores the names of the compositional fields that will
       * be used in the simulation.
       */
      std::vector<std::string> composition_names;

      /**
       * A vector that stores descriptions of each compositional field,
       * including its type (i.e. whether the compositional field corresponds
       * to chemical composition, porosity etc.).
       */
      std::vector<CompositionalFieldDescription> composition_descriptions;

      /**
       * A vector of vectors of composition names that stores the
       * names of the compositional fields corresponding to each field type
       * given in CompositionalFieldDescription.
       */
      std::vector<std::vector<std::string>> composition_names_for_type;

      /**
       * A vector of vectors of composition indices that stores the
       * indices of the compositional fields corresponding to each field type
       * given in CompositionalFieldDescription.
       */
      std::vector<std::vector<unsigned int>> composition_indices_for_type;

      /**
       * List of base element indices used by compositional fields. Cached
       * result returned by get_composition_base_element_indices().
       */
      std::vector<unsigned int> composition_base_element_indices;

      /**
       * Map base_element_index to list of compositional field indices that use
       * that base element. Cached result returned by
       * get_compositional_field_indices_with_base_element();
       */
      std::map<unsigned int, std::vector<unsigned int>> compositional_field_indices_with_base_element;
  };
}


#endif
