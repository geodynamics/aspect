/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/base/std_cxx1x/tuple.h>

namespace aspect
{
  namespace internal
  {

    /**
     * Return pair with @p n_components and a filled ComponentIndices structure.
     */
    template <int dim>
    std::pair<unsigned int, typename Introspection<dim>::ComponentIndices>
    setup_component_indices (const unsigned int n_compositional_fields)
    {
      typename Introspection<dim>::ComponentIndices ci;
      unsigned int comp = 0;

      for (unsigned int i=0; i<dim; ++i)
        ci.velocities[i] = comp++;
      ci.pressure = comp++;
      ci.temperature = comp++;
      for (unsigned int i=0; i<n_compositional_fields; ++i)
        ci.compositional_fields.push_back(comp++);
      return std::make_pair(comp, ci);
    }

    /**
     * Return pair with @p n_blocks and a filled BlockIndices structure.
     */
    template <int dim>
    std::pair<unsigned int, typename Introspection<dim>::BlockIndices>
    setup_blocks (const unsigned int n_compositional_fields,
                  const bool use_direct_solver)
    {
      typename Introspection<dim>::BlockIndices b;

      unsigned int split = (use_direct_solver)?0:1;
      unsigned int block = 0;

      b.velocities = block;
      block += split;
      b.pressure = block++;

      b.temperature = block++;
      for (unsigned int i=0; i<n_compositional_fields; ++i)
        b.compositional_fields.push_back(block++);

      return std::make_pair(block, b);
    }

    /**
     * Return base element structure, FiniteElement spaces, and multiplicities.
     */
    template <int dim>
    std_cxx1x::tuple<typename Introspection<dim>::BaseElements, std::vector<const FiniteElement<dim> *>, std::vector<unsigned int> >
    setup_fes (const Parameters<dim> &parameters)
    {
      typename Introspection<dim>::BaseElements bes;
      std::vector<const FiniteElement<dim> *> fes;
      std::vector<unsigned int> multiplicities;

      unsigned int base_element = 0;

      // u
      fes.push_back(new FE_Q<dim>(parameters.stokes_velocity_degree));
      multiplicities.push_back(dim);
      bes.velocities = base_element++;

      // p
      if (parameters.use_locally_conservative_discretization)
        fes.push_back(new FE_DGP<dim>(parameters.stokes_velocity_degree-1));
      else
        fes.push_back(new FE_Q<dim>(parameters.stokes_velocity_degree-1));
      multiplicities.push_back(1);
      bes.pressure = base_element++;

      // T
      fes.push_back(new FE_Q<dim>(parameters.temperature_degree));
      multiplicities.push_back(1);
      bes.temperature = base_element++;

      // compositions:
      fes.push_back(new FE_Q<dim>(parameters.composition_degree));
      multiplicities.push_back(parameters.n_compositional_fields);
      bes.compositional_fields = base_element++;

      Assert(base_element == fes.size(), ExcInternalError());
      Assert(base_element == multiplicities.size(), ExcInternalError());

      return std_cxx11::make_tuple(bes, fes, multiplicities);
    }

    /**
     * Construct mapping from component to block indices.
     */
    template <int dim>
    std::vector<unsigned int>
    setup_component_to_blocks (const typename Introspection<dim>::ComponentIndices &component_indices,
                               const typename Introspection<dim>::BlockIndices &block_indices,
                               const unsigned int n_components)
    {
      std::vector<unsigned int> components_to_blocks;
      const unsigned int n_compositional_fields = component_indices.compositional_fields.size();

      components_to_blocks.resize(n_components, dealii::numbers::invalid_unsigned_int);
      for (unsigned int d=0; d<dim; ++d)
        components_to_blocks[component_indices.velocities[d]] = block_indices.velocities;
      components_to_blocks[component_indices.pressure] = block_indices.pressure;
      components_to_blocks[component_indices.temperature] = block_indices.temperature;
      for (unsigned int c=0; c<n_compositional_fields; ++c)
        components_to_blocks[component_indices.compositional_fields[c]] = block_indices.compositional_fields[c];

#ifdef DEBUG
      // check we assigned all components
      for (unsigned int c=0; c<n_components; ++c)
        Assert(components_to_blocks[c]!=dealii::numbers::invalid_unsigned_int, ExcInternalError());
#endif

      return components_to_blocks;
    }
  }



  template <int dim>
  Introspection<dim>::Introspection(const Parameters<dim> &parameters)
    :
    n_components (internal::setup_component_indices<dim>(parameters.names_of_compositional_fields.size()).first),
    component_indices (internal::setup_component_indices<dim>(parameters.names_of_compositional_fields.size()).second),
    n_blocks (internal::setup_blocks<dim>(parameters.names_of_compositional_fields.size(), parameters.use_direct_stokes_solver).first),
    block_indices (internal::setup_blocks<dim>(parameters.names_of_compositional_fields.size(), parameters.use_direct_stokes_solver).second),
    extractors (component_indices),
    base_elements (std_cxx1x::get<0>(internal::setup_fes<dim>(parameters))),
    components_to_blocks (internal::setup_component_to_blocks<dim>(component_indices, block_indices, n_components)),
    system_dofs_per_block (n_blocks),
    composition_names(parameters.names_of_compositional_fields),
    fes (std_cxx1x::get<1>(internal::setup_fes<dim>(parameters))),
    multiplicities (std_cxx1x::get<2>(internal::setup_fes<dim>(parameters)))
  {
  }


  template <int dim>
  Introspection<dim>::~Introspection ()
  {
    free_finite_element_data();
  }


  template <int dim>
  void Introspection<dim>::free_finite_element_data ()
  {
    for (unsigned int i=0; i<fes.size(); ++i)
      delete fes[i];
    fes.clear();
    multiplicities.clear();
  }

  namespace
  {
    std::vector<FEValuesExtractors::Scalar>
    make_extractor_sequence (const std::vector<unsigned int> &compositional_fields)
    {
      std::vector<FEValuesExtractors::Scalar> x;
      for (unsigned int i=0; i<compositional_fields.size(); ++i)
        x.push_back (FEValuesExtractors::Scalar(compositional_fields[i]));
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

  template <int dim>
  const std::vector<const FiniteElement<dim> *> &
  Introspection<dim>::get_fes() const
  {
    Assert(fes.size()>0,
           ExcMessage("Error: finite element spaces are only available during construction."));
    return fes;
  }

  template <int dim>
  const std::vector<unsigned int> &
  Introspection<dim>::get_multiplicities() const
  {
    Assert(multiplicities.size()>0,
           ExcMessage("Error: finite element multiplicities are only available during construction."));
    return multiplicities;
  }

}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct Introspection<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
