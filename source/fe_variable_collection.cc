/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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



#include <aspect/global.h>
#include <aspect/fe_variable_collection.h>

namespace aspect
{



  template <int dim>
  VariableDeclaration<dim>::VariableDeclaration(const std::string &name,
                                                const std_cxx11::shared_ptr<FiniteElement<dim> > &fe,
                                                const unsigned int multiplicity,
                                                const unsigned int n_blocks)
    : name(name),
      fe(fe),
      multiplicity(multiplicity),
      n_blocks(n_blocks)
  {
    // TODO: non-scalar FEs are not tested and can not be split into >1 block in general
    Assert(fe->n_components() == 1, ExcNotImplemented());

    Assert(n_blocks == 0
           || n_blocks == 1
           || n_blocks == n_components(),
           ExcMessage("A Variable can only have 0, 1, or n_components() number of blocks."));
  }

  template <int dim>
  VariableDeclaration<dim>::VariableDeclaration()
  {}

  template <int dim>
  VariableDeclaration<dim>::VariableDeclaration(const VariableDeclaration &other)
    : name (other.name),
      fe (other.fe),
      multiplicity (other.multiplicity),
      n_blocks (other.n_blocks)
  {}

  template <int dim>
  unsigned int
  VariableDeclaration<dim>::n_components() const
  {
    return fe->n_components() * multiplicity;
  }


  template <int dim>
  FEVariable<dim>::FEVariable(const VariableDeclaration<dim> &fe_variable,
                              const unsigned int component_index,
                              const unsigned int block_index,
                              const unsigned int base_index)
    : VariableDeclaration<dim> (fe_variable),
      first_component_index (component_index),
      block_index (block_index),
      base_index (base_index),
      scalar_extractor ( (this->n_components()==1) ? component_index : -1),
      vector_extractor ( (this->n_components()==dim) ? component_index : -1)
  {}

  template <int dim>
  const FEValuesExtractors::Scalar &
  FEVariable<dim>::extractor_scalar() const
  {
    Assert(this->n_components()==1,
           ExcMessage("You cannot ask for the scalar extractor of a non-scalar variable."));
    return scalar_extractor;
  }

  template <int dim>
  const FEValuesExtractors::Vector &
  FEVariable<dim>::extractor_vector() const
  {
    Assert(this->n_components()==dim,
           ExcMessage("You cannot ask for the vector extractor of a variable that is not a vector."));
    return vector_extractor;
  }



  template <int dim>
  FEVariableCollection<dim>::FEVariableCollection()
  {}


  template <int dim>
  FEVariableCollection<dim>::FEVariableCollection(const std::vector<VariableDeclaration<dim> > &variable_definitions)
  {
    initialize(variable_definitions);
  }


  template <int dim>
  void
  FEVariableCollection<dim>::initialize(const std::vector<VariableDeclaration<dim> > &variable_definitions)
  {
    variables.clear();
    variables.reserve(variable_definitions.size());

    // store names temporarily to make sure they are unique
    std::set<std::string> names;

    unsigned int component_index = 0;
    unsigned int block_index = 0;

    for (unsigned int i=0; i<variable_definitions.size(); ++i)
      {
        variables.push_back(FEVariable<dim>(variable_definitions[i],
                                            component_index, block_index, i));
        component_index+= variables[i].n_components();
        block_index += variables[i].n_blocks;

        Assert(names.find(variables[i].name)==names.end(),
               ExcMessage("Can not add two variables with the same name."));
        names.insert(variables[i].name);
      }

    Assert(variables.back().n_blocks != 0
           || variables.back().n_components() == 0,
           ExcMessage("The last variable needs to have >0 blocks."));

    n_components_ = component_index;
    n_blocks_ = block_index;

    fes.resize(variables.size());
    multiplicities.resize(variables.size());
    for (unsigned int i=0; i<variables.size(); ++i)
      {
        fes[i] = variables[i].fe.get();
        multiplicities[i] = variables[i].multiplicity;
        variables[i].component_mask = ComponentMask(n_components_, false);
        for (unsigned int c=0; c<variables[i].n_components(); ++c)
          variables[i].component_mask.set(c+variables[i].first_component_index, true);
      }

    components_to_blocks.clear();
    for (unsigned int i=0; i<variables.size(); ++i)
      {
        for (unsigned int c=0; c<variables[i].n_components(); ++c)
          components_to_blocks.push_back(variables[i].block_index
                                         + ((variables[i].n_blocks>1)?c:0));
      }
    Assert(components_to_blocks.size() == n_components(), ExcInternalError());

    Assert(fes.size() == variables.size(), ExcInternalError());
    Assert(multiplicities.size() == variables.size(), ExcInternalError());
  }



  template <int dim>
  const FEVariable<dim> &
  FEVariableCollection<dim>::variable(const std::string &name) const
  {
    for (unsigned int i=0; i<variables.size(); ++i)
      if (variables[i].name == name)
        return variables[i];

    Assert(false, ExcMessage("Variable '" + name + "' not found!"));
    return *(variables.end()); // reference invalid iterator here
  }



  template <int dim>
  bool
  FEVariableCollection<dim>::variable_exists(const std::string &name) const
  {
    for (unsigned int i=0; i<variables.size(); ++i)
      if (variables[i].name == name)
        return true;

    return false;
  }



  template <int dim>
  const std::vector<FEVariable<dim> > &
  FEVariableCollection<dim>::get_variables() const
  {
    return variables;
  }



  template <int dim>
  unsigned int
  FEVariableCollection<dim>::n_components() const
  {
    return n_components_;
  }



  template <int dim>
  unsigned int
  FEVariableCollection<dim>::n_blocks() const
  {
    return n_blocks_;
  }



  template <int dim>
  const std::vector<const FiniteElement<dim> *> &
  FEVariableCollection<dim>::get_fes() const
  {
    return fes;
  }


  template <int dim>
  const std::vector<unsigned int> &
  FEVariableCollection<dim>::get_multiplicities() const
  {
    return multiplicities;
  }


  template <int dim>
  const std::vector<unsigned int> &
  FEVariableCollection<dim>::get_components_to_blocks() const
  {
    return components_to_blocks;
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template struct VariableDeclaration<dim>; \
  template struct FEVariable<dim>; \
  template class FEVariableCollection<dim>; \
   
  ASPECT_INSTANTIATE(INSTANTIATE)
}
