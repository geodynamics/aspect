/*
 Copyright (C) 2023 by the authors of the ASPECT code.
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

#ifndef _aspect_particle_property_elastic_tensor_decomposion_h
#define _aspect_particle_property_elastic_tensor_decomposion_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>
#include <array>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      /**
       * Computes several properties of a elastic tensor stored on a particle.
       * These include the eigenvectors and several projectsion on symmetry axis.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class ElasticTensorDecomposition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * constructor
           */
          ElasticTensorDecomposition();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          void
          initialize () override;

          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * Update function. This function is called every time an update is
           * request by need_update() for every particle for every property.
           *
           * @param [in] data_position An unsigned integer that denotes which
           * component of the particle property vector is associated with the
           * current property. For properties that own several components it
           * denotes the first component of this property, all other components
           * fill consecutive entries in the @p particle_properties vector.
           *
           * @param [in] position The current particle position.
           *
           * @param [in] solution The values of the solution variables at the
           * current particle position.
           *
           * @param [in] gradients The gradients of the solution variables at
           * the current particle position.
           *
           * @param [in,out] particle_properties The properties of the particle
           * that is updated within the call of this function.
           */
          void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim> &position,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim>> &gradients,
                                        const ArrayView<double> &particle_properties) const  override;

          /**
           * This implementation tells the particle manager that
           * we need to update particle properties every time step.
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * Return which data has to be provided to update the property.
           * The integrated strains needs the gradients of the velocity.
           */
          UpdateFlags
          get_needed_update_flags () const  override;

          /**
           * Return a specific permutation based on an index
           */
          static
          std::array<unsigned short int, 3>
          indexed_even_permutation(const unsigned short int index);

          /**
           * Computes the voigt stiffness tensor from the elastic tensor
           */
          static
          SymmetricTensor<2,3>
          compute_voigt_stiffness_tensor(const SymmetricTensor<2,6> &elastic_tensor);

          /**
           * Computes the dilatation stiffness tensor from the elastic tensor
           */
          static
          SymmetricTensor<2,3>
          compute_dilatation_stiffness_tensor(const SymmetricTensor<2,6> &elastic_tensor);

          /**
           * Computes the Symmetry Cartesian Coordinate System (SCCS).
           * With the three SCCS directions, the elastic tensor can be decomposed into the different
           * symmetries in those three SCCS direction, that is, triclinic, monoclinic, orthorhombic,
           * tetragonal, hexagonal, and isotropic (Browaeys & Chevrot, 2004).
           */
          static
          Tensor<2,3> compute_unpermutated_SCCS(const SymmetricTensor<2,3> &dilatation_stiffness_tensor,
                                                const SymmetricTensor<2,3> &voigt_stiffness_tensor);

          /**
           * Uses the Symmetry Cartesian Coordinate System (SCCS) to try the different permutations to
           * determine what is the best projections.
           * This is based on Browaeys and Chevrot (2004), GJI (doi: 10.1111/j.1365-246X.2004.024115.x),
           * which states at the end of paragraph 3.3 that "... an important property of an orthogonal projection
           * is that the distance between a vector $X$ and its orthogonal projection $X_H = p(X)$ on a given
           * subspace is minimum. These two features ensure that the decomposition is optimal once a 3-D Cartesian
           * coordiante systeem is chosen.". The other property they talk about is that "The space of elastic
           * vectors has a finite dimension [...], i.e. using a differnt norm from eq. (2.3 will change disstances
           * but not the resulting decomposition.".
           */
          static
          std::array<std::array<double,3>,7>
          compute_elastic_tensor_SCCS_decompositions(
            const Tensor<2,3> &unpermutated_SCCS,
            const SymmetricTensor<2,6> &elastic_matrix);

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const  override;

          /**
           * Declare paramters
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Parse parameters
           */
          void
          parse_parameters (ParameterHandler &prm) override;


          static SymmetricTensor<2,21> projection_matrix_tric_to_mono;
          static SymmetricTensor<2,9> projection_matrix_mono_to_ortho;
          static SymmetricTensor<2,9> projection_matrix_ortho_to_tetra;
          static SymmetricTensor<2,9> projection_matrix_tetra_to_hexa;
          static SymmetricTensor<2,9> projection_matrix_hexa_to_iso;

        private:
          unsigned int cpo_data_position;
          unsigned int cpo_elastic_tensor_data_position;
      };
    }
  }
}

#endif
