/*
  Copyright (C) 2019 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_rheology_friction_models_h
#define _aspect_material_model_rheology_friction_models_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>
#include <aspect/utilities.h>

#include <aspect/material_model/rheology/drucker_prager.h>
#include <deal.II/fe/component_mask.h>

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {
      /**
       * Enumeration for selecting which type of friction dependence to use.
       *
       * For the type 'static friction', the user-supplied internal angle of friction is used.
       *
       * For the type 'dynamic friction', the friction angle is rate dependent following
       * Equation 13 from \\cite{van_dinther_seismic_2013}.
       *
       * For the type 'differential dynamic friction', the dynamic angle depends on whether
       * flow is convergent, divergent, or transform/neutral according to the
       * projected divergence of the tangential surface velocity.
       */
      enum FrictionMechanism
      {
        static_friction,
        dynamic_friction,
        differential_dynamic_friction,
        function
      };

      /**
       * Numeric codes used to visualize the tectonic regime selected by the
       * differential dynamic friction mechanism.
       */
      enum TectonicRegime
      {
        convergent_regime = -1,
        transform_or_neutral_regime = 0,
        divergent_regime = 1
      };

      template <int dim>
      class FrictionModels : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * A function that computes the new friction angle when it is not independent.
           * Given a volume fraction with the index and a vector of all compositional
           * fields, it returns the newly calculated friction angle.
           */
          double
          compute_friction_angle(const double current_edot_ii,
                                 const unsigned int volume_fraction_index,
                                 const double static_friction_angle,
                                 const Point<dim> &position) const;

          /**
           * A function that returns the selected type of friction dependence.
           */
          FrictionMechanism
          get_friction_mechanism () const;

          /**
           * Return the signed divergence indicator used to select the tectonic
           * regime. This is the previous completed Stokes solution's
           * surface-velocity divergence, projected downward.
           */
          double
          compute_tectonic_divergence_indicator(const Point<dim> &position) const;

          /**
           * Classify the projected surface deformation using the convergence
           * and divergence thresholds.
           */
          TectonicRegime
          compute_tectonic_regime(const Point<dim> &position) const;

        private:
          using SurfaceIndexPoint = boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
          using SurfaceIndexValue = std::pair<SurfaceIndexPoint, double>;
          using SurfaceIndex = boost::geometry::index::rtree<SurfaceIndexValue,
                boost::geometry::index::quadratic<16>>;

          /** Update the cached surface divergence after a nonlinear solve. */
          void
          update_surface_velocity_divergence();

          /** Convert a model position to a depth-independent surface index point. */
          SurfaceIndexPoint
          make_surface_index_point(const Point<dim> &position) const;

          /** Interpolate the cached surface divergence at a lateral position. */
          double
          interpolate_surface_velocity_divergence(const Point<dim> &position) const;

          /**
           * Select the mechanism to be used for the friction dependence.
           * Possible options: static friction | dynamic friction | differential_dynamic_friction | function |
           */
          FrictionMechanism friction_mechanism;

          /**
           * Dynamic friction input parameters
           */

          /**
           * Dynamic angles of internal friction that are used at high strain rates.
           */
          std::vector<double> dynamic_angles_of_internal_friction;

          /**
           * The characteristic strain rate value at which the angle of friction is taken as
           * the mean of the dynamic and the static angle of friction. When the effective
           * strain rate in a cell is very high the dynamic angle of friction is taken, when
           * it is very low the static angle of internal friction is chosen.
           */
          double dynamic_characteristic_strain_rate;

          /**
           * An exponential factor in the equation for the calculation of the friction angle
           * to make the transition between static and dynamic friction angle more smooth or
           * more step-like.
           */
          double dynamic_friction_smoothness_exponent;

          /**
           * Parsed functions that specify the friction angle which must be
           * given in the input file using the function method.
           */
          std::unique_ptr<Functions::ParsedFunction<dim>> friction_function;

          /**
           * The coordinate representation to evaluate the function for the friction angle.
           * Possible choices are depth, cartesian and spherical.
           */
          Utilities::Coordinates::CoordinateSystem coordinate_system_friction_function;

          /**
           * Dynamic angles of internal friction that are used at high strain rates in  converging regions,
           * when using differential dynamic friction.
           */
          std::vector<double> dynamic_angles_of_internal_friction_for_convergence;
          /**
           * Dynamic angles of internal friction that are used at high strain rates in  diverging regions,
           * when using differential dynamic friction.
           */
          std::vector<double> dynamic_angles_of_internal_friction_for_divergence;

          /**
           * Magnitude of negative horizontal divergence required to classify a
           * convergent regime.
           */
          double convergence_threshold;
          /**
           * Positive horizontal divergence required to classify a divergent regime.
           */
          double divergence_threshold;

          /** Maximum depth to which the surface classification is projected. */
          double surface_regime_projection_depth;

          /** Spatial index containing the most recently computed surface divergence. */
          SurfaceIndex surface_divergence_index;
      };
    }
  }
}
#endif
