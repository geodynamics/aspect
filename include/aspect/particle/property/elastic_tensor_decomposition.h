/*
 Copyright (C) 2023 - 2024 by the authors of the ASPECT code.
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

#ifndef _aspect_particle_property_elastic_tensor_decomposition_h
#define _aspect_particle_property_elastic_tensor_decomposition_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      namespace Utilities
      {
        /**
         * Because many places in ASPECT assume that all functions in the namespace
         * <code>aspect::Utilities</code> are available without qualification as
         * <code>Utilities::function</code>, we make sure all these functions
         * are also available inside <code>aspect::Particle::Property::Utilities</code>.
         * This is maybe not the cleanest solution, but it is most compatible
         * with a lot of existing code.
         *
         * We need to do this in every header that creates a new namespace named
         * <code>Utilities</code>, because otherwise the compiler may not find
         * the requested function in the current namespace and issue an error, even
         * though the function is available in the namespace <code>aspect::Utilities</code>.
         */
        using namespace aspect::Utilities;

        /**
         * Return an even permutation based on an index. This is an internal
         * utilities function, also used by the unit tester.
         */
        std::array<unsigned int, 3>
        indexed_even_permutation(const unsigned int index);

        /**
         * Computes the Voigt stiffness tensor from the elastic tensor.
         * The Voigt stiffness tensor (see Browaeys and Chevrot, 2004)
         * defines the stress needed to cause an isotropic strain in the
         * material.
         */
        SymmetricTensor<2,3>
        compute_voigt_stiffness_tensor(const SymmetricTensor<2,6> &elastic_tensor);

        /**
         * Computes the dilatation stiffness tensor from the elastic tensor
         * The dilatational stiffness tensor (see Browaeys and Chevrot, 2004)
         * defines the stress to cause isotropic dilatation in the material.
         */
        SymmetricTensor<2,3>
        compute_dilatation_stiffness_tensor(const SymmetricTensor<2,6> &elastic_tensor);

        /**
         * Computes the Symmetry Cartesian Coordinate System (SCCS).
         *
         * This is based on Browaeys and Chevrot (2004), GJI (doi: 10.1111/j.1365-246X.2004.024115.x),
         * which states at the end of paragraph 3.3 that "... an important property of an orthogonal projection
         * is that the distance between a vector $X$ and its orthogonal projection $X_H = p(X)$ on a given
         * subspace is minimum. These two features ensure that the decomposition is optimal once a 3-D Cartesian
         * coordinate system is chosen.". The other property they talk about is that "The space of elastic
         * vectors has a finite dimension [...], i.e. using a different norm from eq. 2.3 will change distances
         * but not the resulting decomposition.".
         *
         * With the three SCCS directions, the elastic tensor can be decomposed into the different
         * symmetries in those three SCCS direction, that is, triclinic, monoclinic, orthorhombic,
         * tetragonal, hexagonal, and isotropic (Browaeys & Chevrot, 2004).
         *
         * The dilatation_stiffness_tensor defines the stress to cause isotropic dilatation in the material.
         * The voigt_stiffness_tensor defines the stress needed to cause an isotropic strain in the material
         */
        Tensor<2,3> compute_unpermutated_SCCS(const SymmetricTensor<2,3> &dilatation_stiffness_tensor,
                                              const SymmetricTensor<2,3> &voigt_stiffness_tensor);

        /**
         * Uses the Symmetry Cartesian Coordinate System (SCCS) to try the different permutations to
         * determine what is the best projection.
         * This is based on Browaeys and Chevrot (2004), GJI (doi: 10.1111/j.1365-246X.2004.024115.x),
         * which states at the end of paragraph 3.3 that "... an important property of an orthogonal projection
         * is that the distance between a vector $X$ and its orthogonal projection $X_H = p(X)$ on a given
         * subspace is minimum. These two features ensure that the decomposition is optimal once a 3-D Cartesian
         * coordinate system is chosen.". The other property they talk about is that "The space of elastic
         * vectors has a finite dimension [...], i.e. using a different norm from eq. 2.3 will change distances
         * but not the resulting decomposition.".
         */
        std::array<std::array<double,3>,7>
        compute_elastic_tensor_SCCS_decompositions(
          const Tensor<2,3> &unpermutated_SCCS,
          const SymmetricTensor<2,6> &elastic_matrix);


        /**
         * The tensors below can be used to project matrices to different symmetries.
         * See Browaeys and Chevrot, 2004.
         */
        static const SymmetricTensor<2,21> projection_matrix_triclinic_to_monoclinic(
          Tensor<2,21>(
        {
          {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
        })
        );
        static const SymmetricTensor<2,9> projection_matrix_monoclinic_to_orthorhombic(
          Tensor<2,9>(
        {
          {1,0,0,0,0,0,0,0,0},
          {0,1,0,0,0,0,0,0,0},
          {0,0,1,0,0,0,0,0,0},
          {0,0,0,1,0,0,0,0,0},
          {0,0,0,0,1,0,0,0,0},
          {0,0,0,0,0,1,0,0,0},
          {0,0,0,0,0,0,1,0,0},
          {0,0,0,0,0,0,0,1,0},
          {0,0,0,0,0,0,0,0,1}
        })
        );
        static const SymmetricTensor<2,9> projection_matrix_orthorhombic_to_tetragonal(
          Tensor<2,9>(
        {
          {0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
          {0.5,0.5,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
          {0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0},
          {0.0,0.0,0.0,0.5,0.5,0.0,0.0,0.0,0.0},
          {0.0,0.0,0.0,0.5,0.5,0.0,0.0,0.0,0.0},
          {0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0},
          {0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.5,0.0},
          {0.0,0.0,0.0,0.0,0.0,0.0,0.5,0.5,0.0},
          {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0}
        })
        );
        static const SymmetricTensor<2,9> projection_matrix_tetragonal_to_hexagonal(
          Tensor<2,9>(
        {
          {3./8.           , 3./8.           , 0.0, 0.0, 0.0, 1./(4.*sqrt(2.))  , 0.0, 0.0, 1./4.            },
          {3./8.           , 3./8.           , 0.0, 0.0, 0.0, 1./(4.*sqrt(2.))  , 0.0, 0.0, 1./4.            },
          {0.0             , 0.0             , 1.0, 0.0, 0.0, 0.0               , 0.0, 0.0, 0.0              },
          {0.0             , 0.0             , 0.0, 0.5, 0.5, 0.0               , 0.0, 0.0, 0.0              },
          {0.0             , 0.0             , 0.0, 0.5, 0.5, 0.0               , 0.0, 0.0, 0.0              },
          {1./(4.*sqrt(2.)), 1./(4.*sqrt(2.)), 0.0, 0.0, 0.0, 3./4.             , 0.0, 0.0, -1./(2.*sqrt(2.))},
          {0.0             , 0.0             , 0.0, 0.0, 0.0, 0.0               , 0.5, 0.5, 0.0              },
          {0.0             , 0.0             , 0.0, 0.0, 0.0, 0.0               , 0.5, 0.5, 0.0              },
          {1./4.           , 1./4.           , 0.0, 0.0, 0.0, -1./(2.*sqrt(2.)) , 0.0, 0.0, 0.5              }
        })
        );
        static const SymmetricTensor<2,9> projection_matrix_hexagonal_to_isotropic(
          Tensor<2,9>(
        {
          {3./15.      , 3./15.      , 3./15.      , sqrt(2.)/15. , sqrt(2.)/15. , sqrt(2.)/15. , 2./15.       , 2./15.       , 2./15.        },
          {3./15.      , 3./15.      , 3./15.      , sqrt(2.)/15. , sqrt(2.)/15. , sqrt(2.)/15. , 2./15.       , 2./15.       , 2./15.        },
          {3./15.      , 3./15.      , 3./15.      , sqrt(2.)/15. , sqrt(2.)/15. , sqrt(2.)/15. , 2./15.       , 2./15.       , 2./15.        },
          {sqrt(2.)/15., sqrt(2.)/15., sqrt(2.)/15., 4./15.       , 4./15.       , 4./15.       , -sqrt(2.)/15., -sqrt(2.)/15., -sqrt(2.)/15. },
          {sqrt(2.)/15., sqrt(2.)/15., sqrt(2.)/15., 4./15.       , 4./15.       , 4./15.       , -sqrt(2.)/15., -sqrt(2.)/15., -sqrt(2.)/15. },
          {sqrt(2.)/15., sqrt(2.)/15., sqrt(2.)/15., 4./15.       , 4./15.       , 4./15.       , -sqrt(2.)/15., -sqrt(2.)/15., -sqrt(2.)/15. },
          {2./15.      , 2./15.      , 2./15.      , -sqrt(2.)/15., -sqrt(2.)/15., -sqrt(2.)/15., 1./5.        , 1./5.        , 1./5.         },
          {2./15.      , 2./15.      , 2./15.      , -sqrt(2.)/15., -sqrt(2.)/15., -sqrt(2.)/15., 1./5.        , 1./5.        , 1./5.         },
          {2./15.      , 2./15.      , 2./15.      , -sqrt(2.)/15., -sqrt(2.)/15., -sqrt(2.)/15., 1./5.        , 1./5.        , 1./5.         }
        })
        );
      }

      /**
       * Computes several properties of a elastic tensor stored on a particle.
       * These include the eigenvectors and conversions between different forms of 4th order tensors.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class ElasticTensorDecomposition : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor
           */
          ElasticTensorDecomposition() = default;

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
           * @copydoc aspect::Particle::Property::Interface::update_particle_properties()
           */
          void
          update_particle_properties (const ParticleUpdateInputs<dim> &inputs,
                                      typename ParticleHandler<dim>::particle_iterator_range &particles) const override;

          /**
           * This function tells the particle manager that
           * we need to update particle properties.
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

        private:
          unsigned int cpo_elastic_tensor_data_position;
      };
    }
  }
}

#endif
