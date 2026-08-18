/*
 Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_particle_interpolator_quaternion_average_h
#define _aspect_particle_interpolator_quaternion_average_h

#include <aspect/particle/interpolator/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      /**
       * Return the averaged properties of all particles on the given cell.
       *
       * @ingroup ParticleInterpolators
       */
      template <int dim>
      class QuaternionAverage : public Interface<dim>, public aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Return the quaternion average per mineral of all quaternion component fields
           * and return the cell-wise averaged properties of all particles of the cell containing the
           * given positions.
           */
          std::vector<std::vector<double>>
          properties_at_points(const ParticleHandler<dim> &particle_handler,
                               const std::vector<Point<dim>> &positions,
                               const ComponentMask &selected_properties,
                               const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const override;

          // avoid -Woverloaded-virtual:
          using Interface<dim>::properties_at_points;

          /**
           * @copydoc Interface<dim>::declare_parameters()
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * @copydoc Interface<dim>::parse_parameters()
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:

          /**
           * average weighted rotations based on minimizing the frobenius norm between rotations
           * corresponds to R_{mean} = argmin_R(\sum_i||R_i - R||_F)
           * see for reference Markley 2007 (doi:10.2514/1.28949)
           * returns the output as a quaternion in the halfspace where q[0] > 0
           */
          std::array<double,4>
          markley_average(const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const;

          /**
           * returns the geodesic distance between two rotations each represened by a quaternion
           * a geodesic describes the shortest path in a non-euclidian space
           * this metric is equal to the Riemannian metric in SO(3) and lies in the interval [0,pi)
           * for a comparison of different metrics see (https://doi.org/10.1007%2Fs10851-009-0161-2)
           */
          double
          SO3_geodesic(const std::array<double,4> &quaternion1, const std::array<double,4> &quaternion2) const;

          /**
           * sends a quaternion into the fundamental zone of the reference quaternion
           *
           *
           */
          std::array<double,4>
          to_fundamental_zone(const std::array<double,4> &quaternion, const std::array<double,4> &reference_quaternion) const;

          /**
           * average over the reduced space of rotations SO(3)/G
           * where G is the symmetry group G c SO(3) that contains the rotations under which an orthotropic crystal stays invariant. (e.g. Crystallographic texture and group representation, Chapter 6)
           * The average is taken by projecting all
           */
          std::array<double,4>
          symmetry_average(const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const;

          double
          objective_function(const std::array<double,4> &mean_quat, const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const ;

          /**
           * Number of minerals
           */
          unsigned int n_minerals;

          /**
           * the name of the symmetry group of underlying rotations that are supposed to be averaged
           * this determines which symmetry operators are used to project rotations
           * into the fundamental zone of the respective symmetry group
           */
          std::string symmetry_group;

          /**
           * stores the data position of the first quaternion component for each mineral
           * used such that one does not have to sift through all particle properties
           * to find the quaternions.
           */
          std::vector<unsigned int> quaternion_data_pos;

          /**
           * A component mask that determines whether a given particle property is not a quaternion.
           * The properties that are not quaternions and which are selected can be input to the base interpolator.
           */
          ComponentMask not_quaternion_properties;

          /**
           * Scheme that is used to interpolate all particle properties but quaternions
           */
          std::unique_ptr<Interface<dim>> base_interpolator;

      };
    }
  }
}

#endif
