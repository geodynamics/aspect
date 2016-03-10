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

#ifndef __aspect__fe_variable_system_h
#define __aspect__fe_variable_system_h

#include <aspect/global.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe.h>

namespace aspect
{
  using namespace dealii;

  /**
   * A structure to describe everything necessary to define a single variable
   * of the finite element system in isolation. It groups the FiniteElement<dim> with other
   * information like a name and the block structure used in linear systems. A
   * std::vector of these structs will be converted into an ordered collection of
   * objects of type FEVariable represented by FEVariableCollection to form the
   * finite element system.
   *
   * Blocks are used to extract or solve for certain variables separately or
   * in block based preconditioners.  As an example, in a Stokes problem one
   * might want velocity and pressure in a single block (for use with a direct
   * solver), or one block for all velocity components and a second block for
   * the pressure (for a Schur complement-based preconditioner), or separate
   * blocks for pressure and each velocity component.
   */
  template <int dim>
  struct VariableDeclaration
  {
    /**
     * Constructor.
     *
     * @param name A user-friendly and unique string representation of the variable.
     *
     * @param fe The FiniteElement class to use. Currently, this needs to be a scalar
     * finite element class (with exactly one component).
     *
     * @param multiplicity Number of copies of the @p fe to create.
     *
     * @param n_blocks Number of blocks this variable represents inside the
     * linear system. A value of 0 will add the next variable (in the std::vector<VariableDeclaration>
     * that is used to construct the FEVariableCollection) to the current
     * block, while 1 will put the next variable into a new block. A value
     * that is equal to the number of components will create a block for each
     * component.
     */
    VariableDeclaration(const std::string &name,
                        const std_cxx11::shared_ptr<FiniteElement<dim> > &fe,
                        const unsigned int multiplicity,
                        const unsigned int n_blocks);

    /**
     * Default constructor.
     */
    VariableDeclaration();

    /**
     * Copy constructor.
     */
    VariableDeclaration(const VariableDeclaration &other);

    /**
     * Return the total number of components of this variable.
     */
    unsigned int n_components() const;

    /**
     * Name of the variable used in FEVariableCollection to identify it.
     */
    std::string name;

    /**
     * The FiniteElement space.
     */
    std_cxx11::shared_ptr<FiniteElement<dim> > fe;

    /**
     * The multiplicity used in FESystem: how many copies of @p fe are there?
     */
    unsigned int multiplicity;

    /**
     * Number of linear algebra blocks occupied by this variable. See
     * constructor for details.
     */
    unsigned int n_blocks;
  };

  /**
   * Struct to represent a variable as part of a FEVariableCollection. Constructed
   * from an instance of VariableDeclaration<dim> but contains additional
   * information that can be queried.
   */
  template<int dim>
  struct FEVariable: public VariableDeclaration<dim>
  {
      /**
       * Initialize this variable as part of a FESystem with the given
       * indices for @p component, @p block, and @p base.
       */
      FEVariable(const VariableDeclaration<dim> &fe_variable,
                 const unsigned int component_index,
                 const unsigned int block_index,
                 const unsigned int base_index);

      /**
       * The first component index of this variable within the
       * FEVariableCollection.
       */
      unsigned int first_component_index;

      /**
       * The index of the block of the linear system this variable is
       * placed in (within the FEVariableCollection).
       */
      unsigned int block_index;
      /**
       * This index represents the how-manyth variable this is in the
       * FEVariableCollection and consequently the FESystem, so the
       *
       */
      unsigned int base_index;

      /**
       * The component mask of this variable.
       */
      ComponentMask component_mask;

      /**
       * Return a scalar FeValuesExtractor if this variable is scalar.
       */
      const FEValuesExtractors::Scalar &extractor_scalar() const;

      /**
       * Return a vector FeValuesExtractor if this variable is a vector.
       */
      const FEValuesExtractors::Vector &extractor_vector() const;

    private:
      /**
       * Stores the object returned by extractor_scalar().
       */
      FEValuesExtractors::Scalar scalar_extractor;

      /**
       * Stores the object returned by extractor_vector().
       */
      FEValuesExtractors::Vector vector_extractor;
  };


  /**
   * A class that represents a collection of variables (of type FEVariable) to
   * form a Finite Element system.  A dealii::FESystem can be constructed
   * from it using get_fes() and get_multiplicities().  Other information
   * contained in this class on top of FESystem are names, and the block
   * structure for linear systems.
   */
  template <int dim>
  class FEVariableCollection
  {
    public:
      /**
       * Construct an empty object.
       */
      FEVariableCollection();

      /**
       * Construct object from a vector of variables (identical to calling
       * initialize()).
       */
      FEVariableCollection(const std::vector<VariableDeclaration<dim> > &variable_definitions);

      /**
       * Fill this object with the given list of @p variables.
       */
      void initialize(const std::vector<VariableDeclaration<dim> > &variable_definitions);

      /**
       * Return the variable with name @p name. Throws an exception if this
       * variable does not exist.
       */
      const FEVariable<dim> &variable(const std::string &name) const;

      /**
       * Returns true if the variable with @p name exists in the list of
       * variables.
       */
      bool variable_exists(const std::string &name) const;

      /**
       * Return the list of all variables.
       */
      const std::vector<FEVariable<dim> > &get_variables() const;

      /**
       * Return the total number of components in the system.
       */
      unsigned int n_components() const;

      /**
       * Return the total number of block of the system.
       */
      unsigned int n_blocks() const;

      /**
       * Return the vector of finite element spaces used for the construction
       * of the FESystem.
       */
      const std::vector<const FiniteElement<dim> *> &get_fes() const;

      /**
       * Return the vector of multiplicities used for the construction of the
       * FESystem.
       */
      const std::vector<unsigned int> &get_multiplicities() const;

      /**
       * Return a variable that describes for each vector component which
       * vector block it corresponds to.
       */
      const std::vector<unsigned int> &get_components_to_blocks() const;

    protected:
      /**
       * std::vector that contains collection of variables.
       */
      std::vector<FEVariable<dim> > variables;

      /**
       * Total number of components of all variables, returned by n_components().
       */
      unsigned int n_components_;
      /**
       * Total number of blocks of all variables, returned by n_blocks().
       */
      unsigned int n_blocks_;
      /**
       * Data to be used in the FESystem constructor, returned by get_fes().
       */
      std::vector<const FiniteElement<dim> *> fes;
      /**
       * Data to be used in the FESystem constructor, returned by get_multiplicities().
       */
      std::vector<unsigned int> multiplicities;
      /**
       * Mapping from component to block, returned by get_components_to_blocks().
       */
      std::vector<unsigned int> components_to_blocks;
  };

} // namespace aspect

#endif
