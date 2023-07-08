/*
 Copyright (C) 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_property_cpo_elastic_tensor_h
#define _aspect_particle_property_cpo_elastic_tensor_h

#include <aspect/particle/property/interface.h>
#include <aspect/particle/property/crystal_preferred_orientation.h>
#include <aspect/simulator_access.h>
#include <array>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      /**
       * Computes the elastic tensor $C_{ijkl} based on the cpo in both olivine
       * and enstatite. It uses a Voigt average.
       *
       * The layout of the data vector per partcle is the following (note that
       * for this plugin the following dim's are always 3):
       * 1 unrolled tensor -> 3x3x3x3 (dim*dim*dim*dim) doubles, starts at:
       *                                   data_position
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class CpoElasticTensor : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor
           */
          CpoElasticTensor();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          virtual
          void
          initialize ();


          /**
           * This implements the Voigt averaging as described in the equation at the
           * bottom of page 385 in Mainprice (1990):
           * $C^V_{ijkl} = \sum^S_l F_s \sum^{N_s}_l C_{ijkl}/N_s$, where $F_s$ is the
           * grain size, $N_s$ is the number of grains and $C_{ijkl}$ is the elastic
           * tensor. This elastic tensor is computed by the equation above in
           * Mainprice (1990): $C_{ijkl} = R_{ip} R_{jg} R_{kr} R_{is} C_{pgrs}$, where
           * R_{ij} is the cpo orientation matrix.
           */
          virtual
          SymmetricTensor<2,6>
          voigt_average_elastic_tensor (const Particle::Property::CrystalPreferredOrientation<dim> &cpo_particle_property,
                                        const unsigned int cpo_data_position,
                                        const ArrayView<double> &data) const;

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
          virtual
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const;

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
          virtual
          void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim> &position,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim>> &gradients,
                                        const ArrayView<double> &particle_properties) const;

          /**
           * This implementation tells the particle manager that
           * we need to update particle properties every time step.
           */
          UpdateTimeFlags
          need_update () const;

          /**
           * Return which data has to be provided to update the property.
           * The integrated strains needs the gradients of the velocity.
           */
          virtual
          UpdateFlags
          get_needed_update_flags () const;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          virtual
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const;

          /**
           * Loads particle data into variables
           */
          static
          SymmetricTensor<2,6>
          get_elastic_tensor(unsigned int cpo_index,
                             const ArrayView<double> &data);

          /**
           * Stores information in variables into the data array
           */
          static
          void
          set_elastic_tensor(unsigned int cpo_data_position,
                             const ArrayView<double> &data,
                             SymmetricTensor<2,6> &elastic_tensor);

          /**
           * Rotate a 3D 4th order tensor with an other 3D 4th
           */
          static
          Tensor<4,3> rotate_4th_order_tensor(const Tensor<4,3> &input_tensor, const Tensor<2,3> &rotation_tensor);


          /**
           * Rotate a 6x6 voigt matrix with an other 3D 4th
           */
          static
          SymmetricTensor<2,6> rotate_6x6_matrix(const Tensor<2,6> &input_tensor, const Tensor<2,3> &rotation_tensor);

          /**
           * Transform a 4th order tensor into a 6x6 matrix
           */
          static
          SymmetricTensor<2,6> transform_4th_order_tensor_to_6x6_matrix(const Tensor<4,3> &input_tensor);


          /**
           * Transform a 6x6 matrix into a 4th order tensor
           */
          static
          Tensor<4,3> transform_6x6_matrix_to_4th_order_tensor(const SymmetricTensor<2,6> &input_tensor);


          /**
           * Form a 21D vector from a 6x6 matrix
           */
          static
          Tensor<1,21> transform_6x6_matrix_to_21D_vector(const SymmetricTensor<2,6> &input_tensor);

          /**
           * From a 21D vector from a 6xt matrix
           */
          static
          SymmetricTensor<2,6> transform_21D_vector_to_6x6_matrix(const Tensor<1,21> &input_tensor);

          /**
           * Tranform a 4th order tensor directly into a 21D vector.
           */
          static
          Tensor<1,21> transform_4th_order_tensor_to_21D_vector(const Tensor<4,3> &input);

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          unsigned int cpo_data_position;

          SymmetricTensor<2,6> stiffness_matrix_olivine;
          SymmetricTensor<2,6> stiffness_matrix_enstatite;

          unsigned int n_grains;
          unsigned int n_minerals;

          /**
           * The tensor representation of the permutation symbol.
           */
          Tensor<3,3> permutation_operator_3d;

      };
    }
  }
}

#endif
