/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__introspection_h
#define __aspect__introspection_h

#include <deal.II/base/index_set.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe.h>

#include <aspect/fe_variable_collection.h>
#include <aspect/parameters.h>


namespace aspect
{
  using namespace dealii;

  /**
   * Helper function to construct the default list of variables to use
   * based on the given set of @p parameters.
   */
  template <int dim>
  std::vector<VariableDeclaration<dim> >
  construct_default_variables (const Parameters<dim> &parameters);



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
      Introspection (const std::vector<VariableDeclaration<dim> > &variables,
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
       * A variable that holds whether the temperature field should use a
       * discontinuous discretization.
       */
      const bool use_discontinuous_temperature_discretization;

      /**
       * A variable that holds whether the composition field(s) should use a
       * discontinuous discretization.
       */
      const bool use_discontinuous_composition_discretization;

      /**
       * A structure that enumerates the vector components of the finite
       * element that correspond to each of the variables in this problem.
       */
      struct ComponentIndices
      {
        unsigned int       velocities[dim];
        unsigned int       pressure;
        unsigned int       temperature;
        std::vector<unsigned int> compositional_fields;
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
       * If there are compositional fields, they are all discretized with the
       * same base element and, consequently, we only need a single index. If
       * a variable does not exist in the problem (e.g., we do not have
       * compositional fields), then the corresponding index is set to an
       * invalid number.
       */
      struct BaseElements
      {
        unsigned int       velocities;
        unsigned int       pressure;
        unsigned int       temperature;
        unsigned int       compositional_fields;
      };
      /**
       * A variable that enumerates the base elements of the finite element
       * that correspond to each of the variables in this problem.
       */
      const BaseElements base_elements;


      /**
       * A structure that contains component masks for each of the variables
       * in this problem. Component masks are a deal.II concept, see the
       * deal.II glossary.
       */
      struct ComponentMasks
      {
        ComponentMasks (FEVariableCollection<dim> &fevs);

        ComponentMask              velocities;
        ComponentMask              pressure;
        ComponentMask              temperature;
        std::vector<ComponentMask> compositional_fields;
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
         * if velocity and pressure end up in the same block.
         */
        IndexSet locally_owned_pressure_dofs;
      };
      /**
       * A variable that contains index sets describing which of the globally
       * enumerated degrees of freedom are owned by the current processor in a
       * parallel computation.
       */
      IndexSets index_sets;

      /**
       * @}
       */

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
       * A function that gets the name of a compositional field as an input
       * parameter and returns if the compositional field is used in this
       * simulation.
       *
       * @param name The name of compositional field (as specified in the
       * input file)
       */
      bool
      compositional_name_exists (const std::string &name) const;

    private:
      /**
       * A vector that stores the names of the compositional fields that will
       * be used in the simulation.
       */
      std::vector<std::string> composition_names;

  };
}


#endif
