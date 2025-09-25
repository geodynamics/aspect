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


#include <aspect/gravity_model/radial_with_tidal_potential.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/tensor.h>

#include <aspect/gravity_model/radial.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    RadialWithTidalPotential<dim>::gravity_vector (const Point<dim> &p) const
    {
      /**
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

      if (p.norm() == 0.0)
        return Tensor<1,dim>();

      const double r = p.norm();
      return -radialconstant.magnitude * p/r + tidal_gravity;
    }


    template <int dim>
    void
    RadialWithTidalPotential<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial with tidal potential");
        {
          prm.declare_entry ("Mass of perterbing body", "1.898e27",
                             Patterns::Double (),
                             "Mass of body that perturbs gravity of modeled body. "
                             "Default value is for modeling Europa, therefore, mass of Jupiter. "
                             "Units is $kg$.");
          prm.declare_entry ("Semimajor axis of orbit", "670900000",
                             Patterns::Double (),
                             "Length of semimajor axis of orbit between modeled body and perturbing body. "
                             "Default value is for Europa's semimajor axis"
                             "Units is $m$.");
          prm.declare_entry ("Rate of nonsynchronous rotation", "1000",
                             Patterns::Double (),
                             "Rate of nonsynchronous rotation (NSR). "
                             "This works for the modeled body having decoupled rotation between interior layers. "
                             "Default value is for Europa's icy shell. "
                             "This will be converted to angular rate. "
                             "Units is $m/year$");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    RadialWithTidalPotential<dim>::parse_parameters (ParameterHandler &prm)
    {
      AssertThrow (dim==3, ExcMessage ("The 'radial with tidal potential' gravity model "
                                       "can only be used in 3D."));
      
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial with tidal potential");
        {
          M_p = prm.get_double ("Mass of perturbing body");
          a_s = prm.get_double ("Semimajor axis of orbit");
          b_NSR = prm.get_double ("Rate of nonsynchronous rotation");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
      
      radialconstant.initialize_simulator (this->get_simulator());
      radialconstant.parse_parameters(prm);
      radialconstant.initialize();

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
    ASPECT_REGISTER_GRAVITY_MODEL(RadialWithTidalPotential,
                                  "radial with tidal potential",
                                  "A gravity model that is the sum of the `radial constant' model "
                                  "(which is radial, pointing inward if the gravity "
                                  "is positive), "
                                  "and a term that results from a tidal potential and that "
                                  "leads to a gravity field that varies with latitude and longitude."
                                  "The magnitude of gravity for the radial constant part is read from the "
                                  "input file in a section `Gravity model/Radial constant'; the "
                                  "parameters that describe the tidal potential contribution are read "
                                  "from a section ``Gravity model/Radial with tidal potential''.")
  }
}
