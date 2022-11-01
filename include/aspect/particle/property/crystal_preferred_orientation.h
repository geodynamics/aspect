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

#ifndef _aspect_particle_property_cpo_h
#define _aspect_particle_property_cpo_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>
#include <array>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      enum class DeformationType
      {
        passive
      };

      enum class DeformationTypeSelector
      {
        passive
      };

      enum class AdvectionMethod
      {
        forward_euler, backward_euler
      };

      enum class CPODerivativeAlgorithm
      {
        spin_tensor, drex_2004
      };

      /**
       * The plugin manages and computes the evolution of Lattice/Crystal Preferred Orientations (LPO/CPO)
       * on particles. Each ASPECT particle represents many grains. Each grain is assigned a size and a orientation
       * matrix. This allows tracking the LPO evolution with kinematic polycrystal CPO evolution models such
       * as D-Rex (Kaminski and Ribe, 2001; ́Kaminski et al., 2004).
       *
       * This plugin stores M minerals and for each mineral it stores N grains. The total amount of memory stored is 2
       * doubles per mineral plus 10 doubles per grain, resulting in a total memory of M * (2 + N * 10) . The layout of
       * the data stored is the following (note that for this plugin the following dims are always 3):
       * 1. M minerals times
       *    1.1  Mineral deformation type   -> 1 double, at location
       *                                      => 0 + mineral_i * (n_grains * 10 + 2)
       *    2.1. Mineral volume fraction    -> 1 double, at location
       *                                      => 1 + mineral_i * (n_grains * 10 + 2)
       *    2.2. N grains times:
       *         2.1. volume fraction grain -> 1 double, at location:
       *                                      => 2 + grain_i * 10 + mineral_i * (n_grains * 10 + 2)
       *         2.2. rotation_matrix grain -> 9 (Tensor<2,dim>::n_independent_components) doubles, starts at:
       *                                      => 3 + grain_i * 10 + mineral_i * (n_grains * 10 + 2)
       *
       * Last used data entry is n_minerals * (n_grains * 10 + 2).
       *
       * We store the same number of grains for all minerals (e.g. olivine and enstatite
       * grains), although their volume fractions may not be the same. This is because we need a minimum number
       * of grains per tracer to perform reliable statistics on it. This minimum should be the same for all
       * minerals.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class CrystalPreferredOrientation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor
           */
          CrystalPreferredOrientation();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          virtual
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
          virtual
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
          virtual
          void
          update_one_particle_property (const unsigned int data_position,
                                        const Point<dim> &position,
                                        const Vector<double> &solution,
                                        const std::vector<Tensor<1,dim>> &gradients,
                                        const ArrayView<double> &particle_properties) const override;

          /**
           * This implementation tells the particle manager that
           * we need to update particle properties every time step.
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * The CPO of late particles is initialized by interpolating from existing particles.
           */
          InitializationModeForLateParticles
          late_initialization_mode () const override;

          /**
           * Return which data has to be provided to update the property.
           * The integrated strains needs the gradients of the velocity.
           */
          virtual
          UpdateFlags
          get_needed_update_flags () const override;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          virtual
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

          /**
           * Computes the volume fraction and grain orientation derivatives of all the grains of a mineral.
           *
           * @param volume_fractions are the current volume fractions of the grains in a mineral.
           * @param rotation_matrices are the current rotation matrices of the grains in a mineral.
           * @param strain_rate is the strain-rate at the location of the particle.
           * @param velocity_gradient_tensor is the velocity gradient tensor at the location of the particle.
           * @param volume_fraction_mineral is the volume fraction of the current mineral with respect to
           * the other minerals in the particle.
           * @param ref_resolved_shear_stress is the reference resolved shear stress of the mineral.
           */
          std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
          compute_derivatives(const unsigned int cpo_index,
                              const ArrayView<double> &data,
                              const unsigned int mineral_i,
                              const SymmetricTensor<2,3> &strain_rate,
                              const Tensor<2,3> &velocity_gradient_tensor,
                              const std::array<double,4> &ref_resolved_shear_stress) const;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

          /**
           * Return the number of grains per particle
           */
          unsigned int
          get_number_of_grains() const;

          /**
           * Return the number of minerals per particle
           */
          unsigned int
          get_number_of_minerals() const;

          /**
           * @brief Returns the value in the data array representing the deformation type.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to get the value of the deformation type for.
           */
          inline
          double get_deformation_type(const unsigned int cpo_data_position,
                                      const ArrayView<double> &data,
                                      const unsigned int mineral_i) const
          {
            return data[cpo_data_position + 0 + mineral_i * (n_grains * 10 + 2)];
          }

          /**
           * @brief Sets the value in the data array representing the deformation type.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to set the value deformation type for.
           * @param deformation_type The value of the of the deformation type to set.
           */
          inline
          void set_deformation_type(const unsigned int cpo_data_position,
                                    const ArrayView<double> &data,
                                    const unsigned int mineral_i,
                                    const double deformation_type) const
          {
            data[cpo_data_position + 0 + mineral_i * (n_grains * 10 + 2)] = deformation_type;
          }

          /**
           * @brief Returns the value in the data array representing the volume fraction of a mineral.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to get the value of the volume fraction of a mineral for.
           */
          inline
          double get_volume_fraction_mineral(const unsigned int cpo_data_position,
                                             const ArrayView<double> &data,
                                             const unsigned int mineral_i) const
          {
            return data[cpo_data_position + 1 + mineral_i *(n_grains * 10 + 2)];
          }

          /**
           * @brief Sets the value in the data array representing the volume fraction of a mineral.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to set the value of the volume fraction of a mineral for.
           * @param volume_fraction_mineral The value of the of the volume fraction of a mineral to set.
           */
          inline
          void set_volume_fraction_mineral(const unsigned int cpo_data_position,
                                           const ArrayView<double> &data,
                                           const unsigned int mineral_i,
                                           const double volume_fraction_mineral) const
          {
            data[cpo_data_position + 1 + mineral_i *(n_grains * 10 + 2)] = volume_fraction_mineral;
          }

          /**
           * @brief Returns the value in the data array representing the volume fraction of a grain.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to get the value of the volume fraction of a grain for.
           * @param grain_i The grain to get the value of the volume fraction of.
           */
          inline
          double get_volume_fractions_grains(const unsigned int cpo_data_position,
                                             const ArrayView<double> &data,
                                             const unsigned int mineral_i,
                                             const unsigned int grain_i) const
          {
            return data[cpo_data_position + 2 + grain_i * 10 + mineral_i * (n_grains * 10 + 2)];
          }

          /**
           * @brief Sets the value in the data array representing the volume fraction of a grain.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to set the value of the volume fraction of a grain for.
           * @param grain_i The grain to set the value of the volume fraction of.
           * @param volume_fractions_grains The value of the of the volume fraction of a grain to set.
           */
          inline
          void set_volume_fractions_grains(const unsigned int cpo_data_position,
                                           const ArrayView<double> &data,
                                           const unsigned int mineral_i,
                                           const unsigned int grain_i,
                                           const double volume_fractions_grains) const
          {
            data[cpo_data_position + 2 + grain_i * 10 + mineral_i * (n_grains * 10 + 2)] = volume_fractions_grains;
          }

          /**
           * @brief Gets the rotation matrix for a grain in a mineral.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to get the value of the rotation matrix of a grain for.
           * @param grain_i The grain to get the value of the rotation matrix of.
           * @return Tensor<2,3> The rotation matrix of a grain in a mineral
           */
          inline
          Tensor<2,3> get_rotation_matrix_grains(const unsigned int cpo_data_position,
                                                 const ArrayView<double> &data,
                                                 const unsigned int mineral_i,
                                                 const unsigned int grain_i) const
          {
            Tensor<2,3> rotation_matrix;
            for (unsigned int i = 0; i < Tensor<2,3>::n_independent_components ; ++i)
              {
                const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(i);
                rotation_matrix[index] = data[cpo_data_position + 3 + grain_i * 10 + mineral_i * (n_grains * 10 + 2) + i];
              }
            return rotation_matrix;
          }

          /**
           * @brief Sets the rotation matrix for a grain in a mineral.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i The mineral to set the value of the rotation matrix of a grain for.
           * @param grain_i The grain to get the value of the rotation matrix of.
           */
          inline
          void set_rotation_matrix_grains(const unsigned int cpo_data_position,
                                          const ArrayView<double> &data,
                                          const unsigned int mineral_i,
                                          const unsigned int grain_i,
                                          const Tensor<2,3> rotation_matrix) const
          {
            for (unsigned int i = 0; i < Tensor<2,3>::n_independent_components ; ++i)
              {
                const dealii::TableIndices<2> index = Tensor<2,3>::unrolled_to_component_indices(i);
                data[cpo_data_position + 3 + grain_i * 10 + mineral_i * (n_grains * 10 + 2) + i] = rotation_matrix[index];
              }
          }


        private:
          /**
           * Computes a random rotation matrix.
           */
          void
          compute_random_rotation_matrix(Tensor<2,3> &rotation_matrix) const;

          /**
           * @brief Updates the volume fractions and rotation matrices with a Forward Euler scheme.
           *
           * Updates the volume fractions and rotation matrices with a Forward Euler scheme:
           * $x_t = x_{t-1} + dt * x_{t-1} * \frac{dx_t}{dt}$. The function returns the sum of
           * the new volume fractions.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i Which mineral to advect for.
           * @param dt The time step used for the advection step
           * @param derivatives A pair containing the derivatives for the volume fractions and
           * orientations respectively.
           * @return double The sum of all volume fractions.
           */
          double
          advect_forward_euler(const unsigned int cpo_data_position,
                               const ArrayView<double> &data,
                               const unsigned int mineral_i,
                               const double dt,
                               const std::pair<std::vector<double>, std::vector<Tensor<2,3>>> &derivatives) const;

          /**
           * @brief Updates the volume fractions and rotation matrices with a Backward Euler scheme.
           *
           * Updates the volume fractions and rotation matrices with a Backward Euler scheme:
           * $x_t = x_{t-1} + dt * x_{t-1} * \frac{dx_t}{dt}$. The function returns the sum of
           * the new volume fractions.
           *
           * @param cpo_data_position The starting index/position of the cpo data in the particle data vector.
           * @param data The particle data vector.
           * @param mineral_i Which mineral to advect for.
           * @param dt The time step used for the advection step
           * @param derivatives A pair containing the derivatives for the volume fractions and
           * orientations respectively.
           * @return double The sum of all volume fractions.
           */
          double
          advect_backward_euler(const unsigned int cpo_data_position,
                                const ArrayView<double> &data,
                                const unsigned int mineral_i,
                                const double dt,
                                const std::pair<std::vector<double>, std::vector<Tensor<2,3>>> &derivatives) const;


          /**
           * Computes and returns the volume fraction and grain orientation derivatives such that
           * the grains stay the same size and the orientations rotating passively with the particle.
           *
           * @param velocity_gradient_tensor is the velocity gradient tensor at the location of the particle.
           */
          std::pair<std::vector<double>, std::vector<Tensor<2,3>>>
          compute_derivatives_spin_tensor(const Tensor<2,3> &velocity_gradient_tensor) const;

          /**
           * Random number generator used for initalization of particles
           */
          mutable boost::mt19937 random_number_generator;
          unsigned int random_number_seed;

          unsigned int n_grains;

          unsigned int n_minerals;

          /**
           * A vector containing the deformation type selectors provided by the user.
           * Should be one of the following: "Olivine: Karato 2008", "Olivine: A-fabric",
           * "Olivine: B-fabric", "Olivine: C-fabric", "Olivine: D-fabric", "Olivine: E-fabric",
           * "Enstatite" or "Passive".
           */
          std::vector<DeformationTypeSelector> deformation_type_selector;

          /**
           * Store the volume fraction for each mineral.
           */
          std::vector<double> volume_fractions_minerals;

          /**
           * Advection method for particle properties
           */
          AdvectionMethod advection_method;

          /**
           * What algorithm to use to compute the derivatives
           */
          CPODerivativeAlgorithm cpo_derivative_algorithm;

          /**
           * This value determines the tolerance used for the Backward Euler and
           * Crank-Nicolson iterations.
           */
          double property_advection_tolerance;

          /**
           * This value determines the the maximum number of iterations used for the
           * Backward Euler and Crank-Nicolson iterations.
           */
          unsigned int property_advection_max_iterations;

          /**
           * The tensor representation of the permutation symbol.
           */
          Tensor<3,3> permutation_operator_3d;

      };
    }
  }
}

#endif
