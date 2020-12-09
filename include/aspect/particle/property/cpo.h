/*
 Copyright (C) 2020 by the authors of the ASPECT code.

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
        Passive
      };
      enum class DeformationTypeSelector
      {
        Passive
      };

      enum class AdvectionMethod
      {
        ForwardEuler, BackwardEuler, CrankNicolson
      };

      enum class CPODerivativeAlgorithm
      {
        SpinTensor, DRex2004
      };

      /**
       * The plugin manages and computes the evolution of Lattice/Crystal Preferred Orientations (LPO/CPO)
       * on particles. Each ASPECT particle can be assigned many grains. Each grain is assigned a size and a orientation
       * matrix. This allows for LPO evolution tracking with Polycrysal kinematic CPO evolution models such
       * as D-Rex (Kaminski and Ribe, 2001; ́Kaminski et al., 2004).
       *
       * The layout of the data vector per particle is the following (note that for this plugin the following dims are always 3):
       * 1. M mineral times
       *    1.1  Mineral deformation type   -> 1 double, at location
       *                                      => cpo_data_position + 0 + mineral_i * (n_grains * 10 + 2)
       *    2.1. Mineral volume fraction    -> 1 double, at location
       *                                      => cpo_data_position + 1 + mineral_i * (n_grains * 10 + 2)
       *    2.2. N grains times:
       *         2.1. volume fraction grain -> 1 double, at location:
       *                                      => cpo_data_position + 2 + grain_i * 10 + mineral_i * (n_grains * 10 + 2)
       *         2.2. a_cosine_matrix grain -> 9 (Tensor<2,dim>::n_independent_components) doubles, starts at:
       *                                      => cpo_data_position + 3 + grain_i * 10 + mineral_i * (n_grains * 10 + 2)
       *
       * Last used data entry is cpo_data_position  + n_minerals * (n_grains * 10 + 2). In case the Crank-Nicolson
       * advection scheme is used, the previous deriative data is stored behind this. It stores per grain one
       * double for the volume fraction derivative and 9 doubles for the rotation matrix derivative. This results
       * the last data entry in this case being:
       * cpo_data_position + n_minerals * (n_grains * 10 + 2) + mineral_n * (n_grains * 10).
       *
       * We store the same number of grains for all minerals (e.g. olivine and enstatite
       * grains), although their volume fractions may not be the same. This is because we need a minimum number
       * of grains per tracer to perform reliable statistics on it. This minimum should be the same for both
       * olivine and enstatite.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class CPO : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * constructor
           */
          CPO();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run.
           */
          virtual
          void
          initialize ();

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
                                        const std::vector<Tensor<1,dim> > &gradients,
                                        const ArrayView<double> &particle_properties) const;

          /**
           * This implementation tells the particle manager that
           * we need to update particle properties every time step.
           */
          UpdateTimeFlags
          need_update () const;


          InitializationModeForLateParticles
          late_initialization_mode () const;

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
          std::vector<std::pair<std::string, unsigned int> >
          get_property_information() const;

          /**
           * Loads particle data into variables. Intended for use by other parts of aspect.
           */
          static
          void
          load_particle_data(unsigned int cpo_index,
                             const ArrayView<double> &data,
                             std::vector<unsigned int> &deformation_type,
                             std::vector<double> &volume_fraction_mineral,
                             std::vector<std::vector<double>> &volume_fractions_mineral,
                             std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_mineral);


          /**
           * Loads particle data into variables. An extention of load_particle_data which includes the grain
           * derivatvies. Intended mostly for use in this plugin.
           */
          void
          load_particle_data_extended(unsigned int cpo_index,
                                      const ArrayView<double> &data,
                                      std::vector<unsigned int> &deformation_type,
                                      std::vector<double> &volume_fraction_mineral,
                                      std::vector<std::vector<double>> &volume_fractions_mineral,
                                      std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_mineral,
                                      std::vector<std::vector<double> > &volume_fractions_mineral_derivatives,
                                      std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_mineral_derivatives) const;

          /**
           * Stores information in variables into the data array. Intended for use by other parts of ASPECT.
           */
          static
          void
          store_particle_data(unsigned int cpo_data_position,
                              const ArrayView<double> &data,
                              std::vector<unsigned int> &deformation_type,
                              std::vector<double> &volume_fraction_mineral,
                              std::vector<std::vector<double>> &volume_fractions_mineral,
                              std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_mineral);


          /**
           * Stores information in variables into the data array. An extention of load_particle_data which includes the grain
           * derivatvies. Intended mostly for use in this plugin.
           */
          void
          store_particle_data_extended(unsigned int cpo_data_position,
                                       const ArrayView<double> &data,
                                       std::vector<unsigned int> &deformation_type,
                                       std::vector<double> &volume_fraction_mineral,
                                       std::vector<std::vector<double>> &volume_fractions_mineral,
                                       std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_mineral,
                                       std::vector<std::vector<double> > &volume_fractions_mineral_derivatives,
                                       std::vector<std::vector<Tensor<2,3> > > &a_cosine_matrices_mineral_derivatives) const;

          /**
           * Updates the volume fractions and rotation matrices with a Forward Euler scheme.
           */
          double
          advect_forward_euler(std::vector<double> &volume_fractions,
                               std::vector<Tensor<2,3> > &a_cosine_matrices,
                               const std::pair<std::vector<double>, std::vector<Tensor<2,3> > > &derivatives,
                               const double strain_rate_second_invariant,
                               const double dt) const;
          /**
           * Updates the volume fractions and rotation matrices with a Backward Euler scheme.
           */
          double
          advect_backward_euler(std::vector<double> &volume_fractions,
                                std::vector<Tensor<2,3> > &a_cosine_matrices,
                                const std::pair<std::vector<double>, std::vector<Tensor<2,3> > > &derivatives,
                                const double strain_rate_second_invariant,
                                const double dt) const;

          /**
           * Updates the volume fractions and rotation matrices with a Crank-Nicolson scheme.
           */
          double
          advect_Crank_Nicolson(std::vector<double> &volume_fractions,
                                std::vector<Tensor<2,3> > &a_cosine_matrices,
                                const std::pair<std::vector<double>, std::vector<Tensor<2,3> > > &derivatives,
                                std::vector<double> &previous_volume_fraction_derivatives,
                                std::vector<Tensor<2,3> > &previous_a_cosine_matrices_derivatives,
                                const double strain_rate_second_invariant,
                                const double dt) const;


          /**
           * A interface function which delegates the computation the derivatives
           * to specialized functions.
           */
          std::pair<std::vector<double>, std::vector<Tensor<2,3> > >
          compute_derivatives(const std::vector<double> &volume_fractions,
                              const std::vector<Tensor<2,3> > &a_cosine_matrices,
                              const SymmetricTensor<2,3> &strain_rate_nondimensional,
                              const Tensor<2,3> &velocity_gradient_tensor_nondimensional,
                              const double volume_fraction_mineral,
                              const std::array<double,4> &ref_resolved_shear_stress) const;

          /**
           * Computes the update for the property rotating passively with the particle.
           */
          std::pair<std::vector<double>, std::vector<Tensor<2,3> > >
          compute_derivatives_spin_tensor(const std::vector<double> &volume_fractions,
                                          const std::vector<Tensor<2,3> > &a_cosine_matrices,
                                          const SymmetricTensor<2,3> &strain_rate_nondimensional,
                                          const Tensor<2,3> &velocity_gradient_tensor_nondimensional,
                                          const double volume_fraction_mineral,
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
          static
          unsigned int
          get_number_of_grains();

          /**
           * Return the number of minerals per particle
           */
          static
          unsigned int
          get_number_of_minerals();

        private:
          // helpful shortcuts
          double rad_to_degree = 180.0/M_PI;
          double degree_to_rad = M_PI/180.0;

          /**
           * Random number generator used for initalization of particles
           */
          mutable boost::lagged_fibonacci44497 random_number_generator;
          unsigned int random_number_seed;

          static
          unsigned int n_grains;

          static
          unsigned int n_minerals;

          std::vector<DeformationTypeSelector> deformation_type_selector;

          std::vector<double> volume_fractions_minerals;

          /**
           * Advection method for particle properties
           */
          AdvectionMethod advection_method;

          /**
           * What algorithm to use to compute the derivatives
           */
          CPODerivativeAlgorithm cpo_derivative_algorithm;

          double property_advection_tolerance;
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
