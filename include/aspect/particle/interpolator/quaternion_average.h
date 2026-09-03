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
           * Return the average rotation or crystal orientation per mineral of all quaternion component fields.
           * Also returns the interpolated properties of all not quaternion properties using a base interpolation scheme.
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
           * An average for weighted rotations based on minimizing the Frobenius norm between rotations.
           * Corresponds to $\mathbf{R}_{\mathrm{mean}} = \mathrm{argmin}_{\mathbf{R}} (\sum_i ||\mathbf{R}_i - \mathbf{R}||_F)$, represented by a quaternion.
           * Returns the output as a quaternion in the half-space where q.w > 0.
           * See Markley et.al. 2007 (doi:10.2514/1.28949).
           */
          std::array<double,4>
          markley_average(const std::vector<std::array<double,4>> &quat_array, const std::vector<double> &weights) const;

          /**
           * Returns the geodesic distance between two rotations, each represented by a quaternion.
           * A geodesic describes the shortest path in a non-Euclidean space.
           * This metric is equal to the Riemannian metric in $SO(3)$ and lies in the interval $[0,\pi)$.
           * For a comparison of different metrics, see Huynh 2009 (https://doi.org/10.1007%2Fs10851-009-0161-2).
           */
          double
          SO3_geodesic(const std::array<double,4> &quaternion1, const std::array<double,4> &quaternion2) const;

          /**
           * Sends a quaternion into the fundamental zone of the reference quaternion.
           * The fundamental zone contains only one representation of any crystal orientation.
           * In this case, the representation closest to the reference quaternion is chosen.
           * In this function, all possible representations of a crystal orientation for one rotation are hard coded and can be expanded for any other symmetry group.
           */
          std::array<double,4>
          to_fundamental_zone(const std::array<double,4> &quaternion, const std::array<double,4> &reference_quaternion) const;

          /**
           * Average over the reduced space of rotations $SO(3)/G_{cr}$,
           * where $G_{cr} \subset SO(3)$ is the symmetry group that contains the rotations under which an orthotropic crystal stays invariant. (e.g. Crystallographic texture and group representation, Chapter 6)
           * The average is taken by iteratively projecting all rotations into the fundamental zone around a proposed mean value.
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
           * The name of the symmetry group of the rotations that are supposed to be averaged.
           * This determines the symmetry operators that are used to project rotations
           * into a fundamental zone.
           */
          std::string symmetry_group;

          /**
           * Stores the data position of the first quaternion component for each mineral
           * used such that one does not have to sift through all particle properties
           * to find the quaternions.
           */
          std::vector<unsigned int> quaternion_data_pos;

          /**
           * A component mask that determines whether a given particle property is not a quaternion.
           * The properties that are not quaternions and are interpolated from particles to the field are passed to the base interpolator.
           */
          ComponentMask not_quaternion_properties;

          /**
           * Particle interpolater that is used to interpolate all particle properties but quaternions.
           */
          std::unique_ptr<Interface<dim>> base_interpolator;

      };
    }
  }
}

#endif
