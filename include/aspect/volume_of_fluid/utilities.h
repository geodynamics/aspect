/*
 Copyright (C) 2016 - 2024-2018 by the authors of the ASPECT code.

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

#ifndef _aspect_volume_of_fluid_utilities_h
#define _aspect_volume_of_fluid_utilities_h

#include <deal.II/base/point.h>

#include <aspect/global.h>

namespace aspect
{
  namespace VolumeOfFluid
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
       * Function to calculate volume fraction contained by indicator function
       * H(d-normal*(x'-x_{cen}')) on the [0, 1]^dim unit cell where x_{cen} is
       * the unit cell center.
       *
       * Currently only works assuming constant Jacobian determinant.
       */
      double compute_fluid_fraction (const Tensor<1, 2> normal,
                                     const double d);
      double compute_fluid_fraction (const Tensor<1, 3> normal,
                                     const double d);

      /**
       * Function to calculate required value of d to obtain given volume
       * fraction for indicator function H(d-normal*(x'-x_{cen}')) on the [0,
       * 1]^dim unit cell where x_{cen} is the unit cell center.
       *
       * Currently only works assuming constant Jacobian determinant.
       */
      double compute_interface_location (const Tensor<1, 2> normal,
                                         const double volume_fraction);
      double compute_interface_location (const Tensor<1, 3> normal,
                                         const double volume_fraction);

      /**
       * Obtain values at points for a polynomial function that is equivalent to
       * the Heaviside function $H(d-normal*xhat)$ on the unit cell when
       * integrated against polynomials of up to the specified degree.
       *
       * Currently works for degree <=1
       *
       * @param degree Maximum degree for exact integration
       * @param normal Interface normal vector, pointed away from the included region
       * @param d Interface parameter specifying location of interface in unit cell
       * @param points Locations to evaluate the constructed polynomial
       * @param values Values of the constructed polynomial at the specified points
       */
      void xFEM_Heaviside(const unsigned int degree,
                          const Tensor<1, 2> normal,
                          const double d,
                          const std::vector<Point<2>> &points,
                          std::vector<double> &values);
      void xFEM_Heaviside(const unsigned int degree,
                          const Tensor<1, 3> normal,
                          const double d,
                          const std::vector<Point<3>> &points,
                          std::vector<double> &values);

      /**
       * Obtain values at points for a polynomial function that is equivalent to
       * the function $\frac{d}{dd}H(d-normal*xhat)$ on the unit cell when
       * integrated against polynomials of up to the specified degree.
       *
       * Currently works for degree <=1
       *
       * @param degree Maximum degree for exact integration
       * @param normal Interface normal vector, pointed away from the included region
       * @param d Interface parameter specifying location of interface in unit cell
       * @param points Locations to evaluate the constructed polynomial
       * @param values Values of the constructed polynomial at the specified points
       */
      void xFEM_Heaviside_derivative_d(const unsigned int degree,
                                       const Tensor<1, 2> normal,
                                       const double d,
                                       const std::vector<Point<2>> &points,
                                       std::vector<double> &values);
      void xFEM_Heaviside_derivative_d(const unsigned int degree,
                                       const Tensor<1, 3> normal,
                                       const double d,
                                       const std::vector<Point<3>> &points,
                                       std::vector<double> &values);


      /**
       * Function to do Newton iteration calculation of correct d for a given
       * normal to get volume_fraction from xFEM_Heaviside integrated against the given
       * weights.
       *
       * @param degree Maximum degree for exact integration
       * @param normal Interface normal vector, pointed away from the included region
       * @param volume_fraction Cell volume fraction in physical space - used for target value for integration
       * @param vol Cell volume in physical space - used for target value for integration
       * @param epsilon Tolerance for Newton iteration
       * @param points Quadrature points to use for update
       * @param weights JxW values to use for quadrature
       */
      template <int dim>
      double compute_interface_location_newton(const unsigned int degree,
                                               const Tensor<1, dim, double> normal,
                                               const double volume_fraction,
                                               const double vol,
                                               const double epsilon,
                                               const std::vector<Point<dim>> &points,
                                               const std::vector<double> &weights);

      /**
       * Function to calculate volume contained by indicator function
       * $H(d-normal*(x'-x_{cen}'))$ on the $[0, 1]^dim$ unit cell where
       * $x_{cen}$ is the unit cell center, using a polynomial mapping of
       * degree up to "degree".
       *
       * @param degree Maximum degree for exact integration
       * @param normal Interface normal vector, pointed away from the included region
       * @param d Interface parameter specifying location of interface in unit cell by distance from cell center
       * @param points Quadrature points to use for update
       * @param weights JxW values to use for quadrature
       */
      template <int dim>
      double compute_fluid_volume(const unsigned int degree,
                                  const Tensor<1, dim, double> normal,
                                  const double d,
                                  const std::vector<Point<dim>> &points,
                                  const std::vector<double> &weights);

      /**
       * Function to calculate flux volume fraction based on a method of
       * characteristics approximation of the interface on the cell's face over
       * the timestep. Calculation assumes an approximation to the interface of the form
       * $H(d-normal*(x'-x_{face center})-time_direction_derivative*t')$ where
       * $t'$ is in terms of a "unit timestep".
       *
       * @param compute_direction Dimension of unit cell we are currently computing along
       * @param time_direction_derivative Approximated gradient in time for the "level set" describing the interface.
       * @param interface_normal_in_cell The normal vector for the current interface reconstruction in the computing cell.
       * @param d_at_face_center The correct d value to for the interface description on the face we are computing for.
       */
      template <int dim>
      double calculate_volume_flux (const unsigned int compute_direction,
                                    const double time_direction_derivative,
                                    const Tensor<1, dim, double> interface_normal_in_cell,
                                    const double d_at_face_center);
    }
  }
}

#endif
