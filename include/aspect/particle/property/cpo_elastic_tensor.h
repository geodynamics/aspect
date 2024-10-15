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
       * Computes the elastic tensor $C_{ijkl}$ based on the cpo in both olivine
       * and enstatite. It uses a Voigt average.
       *
       * The layout of the data vector per particle is the following (note that
       * for this plugin the following dims are always 3):
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
          void
          initialize () override;

          /**
           * This implements the Voigt averaging as described in the equation at the
           * bottom of page 385 in Mainprice (1990, doi: 10.1016/0098-3004(90)90072-2):
           * $C^V_{ijkl} = \sum^S_l F_s \sum^{N_s}_l C_{ijkl}/N_s$, where $F_s$ is the
           * grain size, $N_s$ is the number of grains and $C_{ijkl}$ is the elastic
           * tensor. This elastic tensor is computed by the equation above in
           * Mainprice (1990): $C_{ijkl} = R_{ip} R_{jg} R_{kr} R_{is} C_{pgrs}$, where
           * R_{ij} is the cpo orientation matrix.
           */
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
           * we need to update particle properties every time step.
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::get_update_flags()
           */
          UpdateFlags
          get_update_flags (const unsigned int component) const override;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

          /**
           * returns the elastic tensor from the particle data
           */
          static
          SymmetricTensor<2,6>
          get_elastic_tensor(unsigned int cpo_index,
                             const ArrayView<double> &data);

          /**
           * Stores the elastic tensor into the particle data array
           */
          static
          void
          set_elastic_tensor(unsigned int cpo_data_position,
                             const ArrayView<double> &data,
                             const SymmetricTensor<2,6> &elastic_tensor);

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Parses the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * The position of the cpo in the particle data array.
           */
          unsigned int cpo_data_position;

          /**
           * The stiffness tensors for olivine and enstatite.
           * Todo: generalize this into a vector.
           */
          SymmetricTensor<2,6> stiffness_matrix_olivine;
          SymmetricTensor<2,6> stiffness_matrix_enstatite;

          /**
           * The number of grains per particle.
           */
          unsigned int n_grains;

          /**
           * The number of minerals per particle.
           */
          unsigned int n_minerals;

      };
    }
  }
}

#endif
