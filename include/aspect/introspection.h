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
/*  $Id$  */


#ifndef __aspect__introspection_h
#define __aspect__introspection_h

#include <deal.II/base/index_set.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_values_extractors.h>


namespace aspect
{
  using namespace dealii;

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
  struct Introspection
  {
    public:
      /**
       * Constructor.
       * @param n_compositional_fields The number of compositional fields that
       * will be used in this simulation. This is used in initializing the
       * fields of this class.
       */
      Introspection (const unsigned int n_compositional_fields);

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
       * The number of vector blocks. This equals $3+n_c$ where, in comparison
       * to the n_components field, the velocity components form a single
       * block.
       */
      const unsigned int n_blocks;

      /**
       * A structure that contains FEValues extractors for every block of the
       * finite element used in the overall description.
       */
      struct Extractors
      {
        Extractors (const unsigned int n_compositional_fields);

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
       * A structure that enumerates the vector components of the finite
       * element that correspond to each of the variables in this problem.
       */
      struct ComponentIndices
      {
        ComponentIndices (const unsigned int n_compositional_fields);

        static const unsigned int       velocities[dim];
        static const unsigned int       pressure    = dim;
        static const unsigned int       temperature = dim+1;
        const std::vector<unsigned int> compositional_fields;
      };
      /**
       * A variable that enumerates the vector components of the finite
       * element that correspond to each of the variables in this problem.
       */
      ComponentIndices component_indices;

      /**
       * A structure that enumerates the vector blocks of the finite element
       * that correspond to each of the variables in this problem.
       */
      struct BlockIndices
      {
        BlockIndices (const unsigned int n_compositional_fields);

        static const unsigned int       velocities  = 0;
        static const unsigned int       pressure    = 1;
        static const unsigned int       temperature = 2;
        const std::vector<unsigned int> compositional_fields;
      };
      /**
       * A variable that enumerates the vector blocks of the finite element
       * that correspond to each of the variables in this problem.
       */
      BlockIndices block_indices;


      /**
       * A structure that contains component masks for each of the variables
       * in this problem. Component masks are a deal.II concept, see the
       * deal.II glossary.
       */
      struct ComponentMasks
      {
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
      ComponentMasks component_masks;

      /**
       * A variable that describes for each vector component which vector
       * block it corresponds to.
       */
      const std::vector<unsigned int> components_to_blocks;

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
  };
}


#endif
