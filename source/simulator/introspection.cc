/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/introspection.h>
#include <aspect/global.h>

namespace aspect
{
  namespace
  {
    template <int dim>
    std::vector<unsigned int>
    component_to_block_mapping (const unsigned int n_components,
                               const bool split_vel_pressure,
                               const bool add_compaction_pressure)
    {
      // set up a mapping between vector components to the blocks they
      // correspond to.
      std::vector<unsigned int> components_to_blocks (n_components, 0U);

      // velocity is always 0, so start at pressure:
      unsigned int start_idx = dim;

      if (!split_vel_pressure)
        ++start_idx; // skip pressure, so it will be block 0

      if (!split_vel_pressure && add_compaction_pressure)
        ++start_idx; // skip compaction pressure, so it will be block 0

      // number the remainder increasing from 1:
      unsigned int block = 0;
      for (unsigned int i=start_idx; i < n_components; ++i)
        components_to_blocks[i] = (++block);

      return components_to_blocks;
    }
  }


  template <int dim>
  Introspection<dim>::Introspection(const bool split_vel_pressure,
      const bool add_compaction_pressure,
                                    const std::vector<std::string> &names_of_compositional_fields)
    :
    n_components (dim+2+names_of_compositional_fields.size()+(add_compaction_pressure?1:0)),
    n_blocks (((split_vel_pressure)?3:2)+names_of_compositional_fields.size()),
    extractors (names_of_compositional_fields.size(), add_compaction_pressure),
    component_indices (names_of_compositional_fields.size(), add_compaction_pressure),
    block_indices (names_of_compositional_fields.size(), split_vel_pressure),
    base_elements (names_of_compositional_fields.size(), add_compaction_pressure),
    components_to_blocks (component_to_block_mapping<dim>(n_components, split_vel_pressure, add_compaction_pressure)),
    system_dofs_per_block (n_blocks),
    composition_names(names_of_compositional_fields)
  {}


  namespace
  {
    std::vector<unsigned int>
    half_open_sequence (const unsigned int begin,
                        const unsigned int end)
    {
      std::vector<unsigned int> x;
      for (unsigned int i=begin; i<end; ++i)
        x.push_back (i);
      return x;
    }
  }


  template <int dim>
  Introspection<dim>::ComponentIndices::
  ComponentIndices (const unsigned int n_compositional_fields,
                    const bool add_compaction_pressure)
    :
    pressure (dim),
    compaction_pressure (add_compaction_pressure ? dim+1 : numbers::invalid_unsigned_int),
    temperature (add_compaction_pressure ? dim+2 : dim+1),
    compositional_fields (half_open_sequence(temperature+1,
                                             temperature+1+n_compositional_fields))
  {
    for (unsigned int i=0;i<dim;++i)
      velocities[i]=i;
  }


  template <int dim>
  Introspection<dim>::BlockIndices::
  BlockIndices (const unsigned int n_compositional_fields,
                const bool split_vel_pressure)
    :
    velocities(0),
    pressure (split_vel_pressure?1:0),
    compaction_pressure (split_vel_pressure?1:0),
    temperature (split_vel_pressure?2:1),
    compositional_fields (half_open_sequence(
                            (split_vel_pressure?3:2),
                            (split_vel_pressure?3:2)+n_compositional_fields))
  {}


  template <int dim>
  Introspection<dim>::BaseElements::
  BaseElements (const unsigned int n_compositional_fields,
                const bool add_compaction_pressure)
    :
    velocities(0),
    pressure (1),
    compaction_pressure (add_compaction_pressure ? 1 : numbers::invalid_unsigned_int),
    temperature (2),
    compositional_fields (n_compositional_fields > 0 ? 3 : numbers::invalid_unsigned_int)
  {}


  namespace
  {
    std::vector<FEValuesExtractors::Scalar>
    half_open_extractor_sequence (const unsigned int begin,
                                  const unsigned int end)
    {
      std::vector<FEValuesExtractors::Scalar> x;
      for (unsigned int i=begin; i<end; ++i)
        x.push_back (FEValuesExtractors::Scalar(i));
      return x;
    }
  }

  template <int dim>
  Introspection<dim>::Extractors::Extractors (const unsigned int n_compositional_fields,
                                              const bool add_compaction_pressure)
    :
    velocities (0),
    pressure (dim),
    compaction_pressure (add_compaction_pressure ? dim+1 : numbers::invalid_unsigned_int),
    temperature (add_compaction_pressure ? dim+2 : dim+1),
    compositional_fields (half_open_extractor_sequence (dim+2+(add_compaction_pressure?1:0),
                                                        dim+2+(add_compaction_pressure?1:0)+n_compositional_fields))
  {}

  template <int dim>
  unsigned int
  Introspection<dim>::compositional_index_for_name (const std::string &name) const
  {
    std::vector<std::string>::const_iterator it = std::find(composition_names.begin(), composition_names.end(), name);
    if (it == composition_names.end())
      {
        AssertThrow (false, ExcMessage ("The compositional field " + name +
                                        " you asked for is not used in the simulation."));
      }
    else
      return it - composition_names.begin();
    return numbers::invalid_unsigned_int;
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
  bool
  Introspection<dim>::compositional_name_exists (const std::string &name) const
  {
    return (std::find(composition_names.begin(), composition_names.end(), name) != composition_names.end()
            ?
            true
            :
            false);
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Introspection<dim>;\
   
  ASPECT_INSTANTIATE(INSTANTIATE)
}
