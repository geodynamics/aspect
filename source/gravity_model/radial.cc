/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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


#include <aspect/gravity_model/radial.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/tensor.h>

#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

namespace aspect
{
  namespace GravityModel
  {
// ------------------------------ RadialConstant -------------------
    template <int dim>
    Tensor<1,dim>
    RadialConstant<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double r = p.norm();
      return -magnitude * p/r;
    }



    template <int dim>
    void
    RadialConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          prm.declare_entry ("Magnitude", "9.81",
                             Patterns::Double (),
                             "Magnitude of the gravity vector in $m/s^2$. For positive values "
                             "the direction is radially inward towards the center of the earth.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial constant");
        {
          magnitude = prm.get_double ("Magnitude");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
                   this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::ellipsoidal,
                   ExcMessage ("Gravity model 'radial constant' should not be used with geometry models that "
                               "do not have either a spherical or ellipsoidal natural coordinate system."));
    }

// ------------------------------ RadialEarthLike -------------------

    template <int dim>
    void
    RadialEarthLike<dim>::initialize ()
    {
      AssertThrow(false,
                  ExcMessage("The 'radial earth-like' gravity model has been removed "
                             "due to its misleading name. The available AsciiData gravity "
                             "model (using default parameters) is much more earth-like, since "
                             "it uses the gravity profile used in the construction of the "
                             "Preliminary Reference Earth Model (PREM, Dziewonski and Anderson, "
                             "1981). Use the 'ascii data' model instead of 'radial earth-like'."));
    }



    template <int dim>
    Tensor<1,dim>
    RadialEarthLike<dim>::gravity_vector (const Point<dim> &/*p*/) const
    {
      // We should never get here, because of the assertion in initialize().
      AssertThrow(false,ExcInternalError());
      return Tensor<1,dim>();
    }


// ----------------------------- RadialLinear ----------------------


    template <int dim>
    Tensor<1,dim>
    RadialLinear<dim>::gravity_vector (const Point<dim> &p) const
    {
      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double depth = this->get_geometry_model().depth(p);
      const double maximal_depth = this->get_geometry_model().maximal_depth();

      return  (-magnitude_at_surface * (1.0 - depth/maximal_depth)
               -magnitude_at_bottom  * (depth/maximal_depth))
              * p/p.norm();
    }


    template <int dim>
    void
    RadialLinear<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear");
        {
          prm.declare_entry ("Magnitude at surface", "9.8",
                             Patterns::Double (),
                             "Magnitude of the radial gravity vector "
                             "at the surface of the domain. "
                             "Units: \\si{\\meter\\per\\second\\squared}.");
          prm.declare_entry ("Magnitude at bottom", "10.7",
                             Patterns::Double (),
                             "Magnitude of the radial gravity vector "
                             "at the bottom of the domain. `Bottom' means the"
                             "maximum depth in the chosen geometry, and for "
                             "example represents the core-mantle boundary in "
                             "the case of the `spherical shell' geometry model, "
                             "and the center in the case of the `sphere' "
                             "geometry model. Units: \\si{\\meter\\per\\second\\squared}.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialLinear<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear");
        {
          magnitude_at_surface = prm.get_double ("Magnitude at surface");
          magnitude_at_bottom = prm.get_double ("Magnitude at bottom");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical ||
                   this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::ellipsoidal,
                   ExcMessage ("Gravity model 'radial linear' should not be used with geometry models that "
                               "do not have either a spherical or ellipsoidal natural coordinate system."));

    }



// ----------------------------- RadialLinearWithTidalPotential ----------------------


    template <int dim>
    Tensor<1,dim>
    RadialLinearWithTidalPotential<dim>::gravity_vector (const Point<dim> &p) const
    {
      /**
       * This gravity model is a plugin for a case where the
       * tidal potential by flattening and non-synchronnous rotation changes gravity with position and time.
       *
       * The equation implemented in this heating model is from Tobie et al. (2025) (https://doi.org/10.1007/s11214-025-01136-y),
       * which is defined as:
       * g = -magnitude - gradient (-tidal potential)).
       * potential = (3 G M_p) / (2 a_s^3) * r^2 * (Tstar + T0)
       * Tstar = 1/6 *(1-3*cos(theta)^2) and T0=1/2sin(theta)^2*cos(2*lambda + 2*b*t)
       * where G = gravitational constant, M_p = mass of planet, a_s = semimajor axis of satellite's orbit, b = angular rate of nonsynchronous rotation.
       * r, theta and lambda are radial distance, polar angle and azimuthal angle, respectively.
       * b [1/s] = b_NSR [m/yr] / (circumference of satellite [m]) / year_to_seconds [s/yr] * 2 * pi
       *
       * Notation of this potential equation is converted from spherical coordinates to cartesian coordinates.
       * Therefore, gradient of potential is (3 G M_p) / (2 a_s^3) * ( 1 / 6 * ( x^2 + y^2 - 2 * z^2) + 1 / 2 * (C1*(x^2 + y^2) - 2 * C2 * x * y)))
       * where C1 = cos(2*b*t) and C2 = sin(2*b*t)
       */
      const Tensor<1,dim> e_x = Point<dim>::unit_vector(0);
      const Tensor<1,dim> e_y = Point<dim>::unit_vector(1);
      const Tensor<1,dim> e_z = Point<dim>::unit_vector(2);

      const double t = (this->simulator_is_past_initialization()) ? this->get_time() : 0.0;
      const double x = p[0];
      const double y = p[1];
      const double z = p[2];

      const double C1 = std::cos( 2. * b_NSR / R1 / year_in_seconds * t);
      const double C2 = std::sin( 2. * b_NSR / R1 / year_in_seconds * t);

      const double dTstar_over_dx = 1. / 6. * ( 2. * x );
      const double dT0_over_dx    = 1. / 2. * ( 2. * C1 * x - 2. * C2 * y );

      const double dTstar_over_dy = 1. / 6. * ( 2. * y );
      const double dT0_over_dy    = 1. / 2. * ( -2. * C1 * y - 2. * C2 * x );

      const double dTstar_over_dz = 1. / 6. * ( -4. * z );
      const double dT0_over_dz    = 0;

      const double G = aspect::constants::big_g;
      const double T_factor = 3. * G * M_p / 2. / a_s / a_s / a_s;

      const Tensor<1,dim> tidal_gravity = T_factor *
                                          ( (dTstar_over_dx + dT0_over_dx) * e_x
                                            + (dTstar_over_dy + dT0_over_dy) * e_y
                                            + (dTstar_over_dz + dT0_over_dz) * e_z);

      return RadialLinear<dim>::gravity_vector(p) + tidal_gravity;
    }


    template <int dim>
    void
    RadialLinearWithTidalPotential<dim>::declare_parameters (ParameterHandler &prm)
    {
      RadialLinear<dim>::declare_parameters (prm);

      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear with tidal potnetial");
        {
          prm.declare_entry ("Mass of planet", "1.898e27",
                             Patterns::Double (),
                             "Mass of satellte's host planet. "
                             "Default value is for Europa, therefore, mass of Jupiter. "
                             "Units is $kg$.");
          prm.declare_entry ("Semimajor axis of satellite", "670900000",
                             Patterns::Double (),
                             "Length of semimajor axis of satellite. "
                             "Default value is for Europa's semimajor axis"
                             "Units is $m$.");
          prm.declare_entry ("Rate of nonsynchronous rotation", "1000",
                             Patterns::Double (),
                             "Rate of nonsynchronous rotation (NSR). "
                             "This will be converted to angular rate. "
                             "Units is $m/year$");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialLinearWithTidalPotential<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (dim==3, ExcMessage ("The 'radial linear with tidal potential' gravity model "
                                       "can only be used in 3D."));
      this->RadialLinear<dim>::parse_parameters (prm);

      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial linear with tidal potnetial");
        {
          M_p = prm.get_double ("Mass of planet");
          a_s = prm.get_double ("Semimajor axis of satellite");
          b_NSR = prm.get_double ("Rate of nonsynchronous rotation");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
      // This effect of tidal potential only works if the geometry is derived from
      // a spherical model (i.e. a sphere, spherical shell or chunk)
      if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()))
        {
          R1 = Plugins::get_plugin_as_type<const GeometryModel::Sphere<dim>>
               (this->get_geometry_model()).radius();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()))
        {
          R1 = Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>>
               (this->get_geometry_model()).outer_radius();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>>(this->get_geometry_model()))
        {
          R1 = Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>>
               (this->get_geometry_model()).outer_radius();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::EllipsoidalChunk<dim>>(this->get_geometry_model()))
        {
          const auto &gm = Plugins::get_plugin_as_type<const GeometryModel::EllipsoidalChunk<dim>>
                           (this->get_geometry_model());
          // TODO
          // If the eccentricity of the EllipsoidalChunk is non-zero, the radius can vary along a boundary,
          // but the maximal depth is the same everywhere and we could calculate a representative pressure
          // profile. However, it requires some extra logic with ellipsoidal
          // coordinates, so for now we only allow eccentricity zero.
          // Using the EllipsoidalChunk with eccentricity zero can still be useful,
          // because the domain can be non-coordinate parallel.

          AssertThrow(gm.get_eccentricity() == 0.0,
                      ExcNotImplemented("This plugin cannot be used with a non-zero eccentricity. "));

          R1 = gm.get_semi_major_axis_a();
        }
      else
        {
          Assert (false, ExcMessage ("This initial condition can only be used if the geometry "
                                     "is a sphere, a spherical shell, a chunk or an "
                                     "ellipsoidal chunk."));
          R1 = numbers::signaling_nan<double>();
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace GravityModel
  {
    ASPECT_REGISTER_GRAVITY_MODEL(RadialConstant,
                                  "radial constant",
                                  "A gravity model in which the gravity has a constant magnitude "
                                  "and the direction is radial (pointing inward if the value "
                                  "is positive). The magnitude is read from the parameter "
                                  "file in subsection 'Radial constant'.")

    ASPECT_REGISTER_GRAVITY_MODEL(RadialEarthLike,
                                  "radial earth-like",
                                  "This plugin has been removed due to its misleading name. "
                                  "The included profile was hard-coded and was less earth-like "
                                  "than the `ascii data' plugin, which uses the profile "
                                  "of the Preliminary Reference Earth Model (PREM). Use `ascii data' "
                                  "instead of `radial earth-like'.")

    ASPECT_REGISTER_GRAVITY_MODEL(RadialLinear,
                                  "radial linear",
                                  "A gravity model which is radial (pointing inward if the gravity "
                                  "is positive) and the magnitude changes linearly with depth. The "
                                  "magnitude of gravity at the surface and bottom is read from the "
                                  "input file in a section ``Gravity model/Radial linear''.")

    ASPECT_REGISTER_GRAVITY_MODEL(RadialLinearWithTidalPotential,
                                  "radial linear with tidal potential",
                                  "A gravity model that is the sum of the `radial linear' model "
                                  "(which is radial, pointing inward if the gravity "
                                  "is positive, and a magnitude that changes linearly with depth), "
                                  "and a term that results from a tidal gravity potential and that "
                                  "leads to a gravity field that varies with latitude and longitude."
                                  "The magnitude of gravity for the linear radial part is read from the "
                                  "input file in a section `Gravity model/Radial linear'; the "
                                  "parameters that describe the tidal potential contribution are read "
                                  "from a section ``Gravity model/Radial linear with tidal potential''.")
  }
}
