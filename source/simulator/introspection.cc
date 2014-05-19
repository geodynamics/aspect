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
  template <>
  const unsigned int
  Introspection<2>::ComponentIndices::velocities[2] = { 0, 1 };

  template <>
  const unsigned int
  Introspection<3>::ComponentIndices::velocities[3] = { 0, 1, 2 };

  template <int dim>
  const unsigned int
  Introspection<dim>::ComponentIndices::pressure;

  template <int dim>
  const unsigned int
  Introspection<dim>::ComponentIndices::temperature;

  template <int dim>
  const unsigned int
  Introspection<dim>::BaseElements::velocities;

  template <int dim>
  const unsigned int
  Introspection<dim>::BaseElements::pressure;

  template <int dim>
  const unsigned int
  Introspection<dim>::BaseElements::temperature;


  namespace
  {
    template <int dim>
    std::vector<unsigned int>
    component_to_block_mapping (const unsigned int n_components,
        const bool split_vel_pressure)
    {
      // set up a mapping between vector components to the blocks they
      // correspond to.
      std::vector<unsigned int> components_to_blocks (n_components, 0U);
      if (split_vel_pressure)
        {
          // each variable has its own block except
          // for the velocities which are all mapped into block 0
          components_to_blocks[dim] = 1;
          components_to_blocks[dim+1] = 2;
          for (unsigned int i=dim+2; i<n_components; ++i)
            components_to_blocks[i] = i-dim+1;
        }
      else
        {
          // here velocity and pressure is block 0:
          components_to_blocks[dim+1] = 1;
          for (unsigned int i=dim+2; i<n_components; ++i)
            components_to_blocks[i] = i-dim;
        }

      return components_to_blocks;
    }
  }


  template <int dim>
  Introspection<dim>::Introspection(const unsigned int n_compositional_fields,
      const bool split_vel_pressure)
    :
    n_components (dim+2+n_compositional_fields),
    n_blocks (((split_vel_pressure)?3:2)+n_compositional_fields),
    extractors (n_compositional_fields),
    component_indices (n_compositional_fields),
    block_indices (n_compositional_fields, split_vel_pressure),
    base_elements (n_compositional_fields),
    components_to_blocks (component_to_block_mapping<dim>(n_components, split_vel_pressure)),
    system_dofs_per_block (n_blocks)
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
  ComponentIndices (const unsigned int n_compositional_fields)
    :
    compositional_fields (half_open_sequence(dim+2, dim+2+n_compositional_fields))
  {}


  template <int dim>
  Introspection<dim>::BlockIndices::
  BlockIndices (const unsigned int n_compositional_fields,
      const bool split_vel_pressure)
    :
    velocities(0),
    pressure (split_vel_pressure?1:0),
    temperature (split_vel_pressure?2:1),
    compositional_fields (half_open_sequence(
        (split_vel_pressure?3:2),
        (split_vel_pressure?3:2)+n_compositional_fields))
  {}


  template <int dim>
  Introspection<dim>::BaseElements::
  BaseElements (const unsigned int n_compositional_fields)
    :
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
  Introspection<dim>::Extractors::Extractors (const unsigned int n_compositional_fields)
    :
    velocities (0),
    pressure (dim),
    temperature (dim+1),
    compositional_fields (half_open_extractor_sequence (dim+2, dim+2+n_compositional_fields))
  {
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Introspection<dim>;\
   
  ASPECT_INSTANTIATE(INSTANTIATE)
}
