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
   * The introspection class provides information about the simulation
   * as a whole. In particular, it provides symbolic names for things
   * like the velocity, pressure and other variables, along with their
   * corresponding vector and scalar FEValues extractors, component
   * masks, etc.
   *
   * The purpose of this class is primarily to provide these symbolic
   * names so that we do not have to use implicit knowledge about the
   * ordering of variables (e.g., earlier versions of ASPECT had
   * many places where we built a scalar FEValues extractor at
   * component 'dim' since that is where we knew that the pressure
   * lies in the finite element; this kind of implicit knowledge is
   * no longer necessary with the Introspection class). The Introspection
   * class is used both by the Simulator class itself, but is also
   * exported to plugins via the SimulatorAccess class.
   *
   * @ingroup Simulator
   */
  template <int dim>
  struct Introspection
  {
    public:
      Introspection (const unsigned int n_compositional_fields);

      /**
       * @name Things that are independent of the current mesh
       * @{
       */
      const unsigned int n_components;
      const unsigned int n_blocks;

      struct Extractors
      {
        Extractors (const unsigned int n_compositional_fields);

        FEValuesExtractors::Vector              velocities;
        FEValuesExtractors::Scalar              pressure;
        FEValuesExtractors::Scalar              temperature;
        std::vector<FEValuesExtractors::Scalar> compositional_fields;
      };
      Extractors extractors;


      struct ComponentMasks
      {
        ComponentMask              velocities;
        ComponentMask              pressure;
        ComponentMask              temperature;
        std::vector<ComponentMask> compositional_fields;
      };
      ComponentMasks component_masks;

      std::vector<unsigned int> components_to_blocks;

      /**
       * @}
       */

      /**
       * @name Things that depend on the current mesh
       * @{
       */
      std::vector<unsigned int> system_dofs_per_block;

      struct IndexSets
      {
        IndexSet system_relevant_set;

        std::vector<IndexSet> system_partitioning;
        std::vector<IndexSet> system_relevant_partitioning;
      };
      IndexSets index_sets;

      /**
       * @}
       */
  };
}


#endif
