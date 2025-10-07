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

namespace aspect
{
  namespace GravityModel
  {
    template <int dim>
    Tensor<1,dim>
    RadialWithTidalPotential<dim>::gravity_vector (const Point<dim> &p) const
    {
    // This plugin is not implemented for 2D models
    AssertThrow(false, ExcNotImplemented());
    return Tensor<1,dim>();
    }

    template <>
    Tensor<1,3>
    RadialWithTidalPotential<3>::gravity_vector (const Point<3> &p) const
    {
      const unsigned int dim = 3;
      /**
       * Notation of this potential equation is converted from spherical coordinates to cartesian coordinates.
       * Therefore, gradient of potential is (3 G M_p) / (2 a_s^3) * ( 1 / 6 * ( x^2 + y^2 - 2 * z^2) + 1 / 2 * (C1*(x^2 + y^2) - 2 * C2 * x * y)))
       * where C1 = cos(2*b*t) and C2 = sin(2*b*t)
       */
      const double t = (this->simulator_is_past_initialization()) ? this->get_time() : 0.0;

      const double C1 = std::cos( 2. * b * t);
      const double C2 = std::sin( 2. * b * t);

      const Tensor<1,dim> dTstar_gradient ({1./3. * p[0], 1./3. * p[1], -2./3. * p[2]});

      const Tensor<1,dim> dT0_gradient ({C1*p[0] - C2*p[1], -C1*p[1] - C2*p[0], 0});

      const double G = aspect::constants::big_g;
      const double T_factor = 3. * G * M_p / ( 2. * a_s * a_s * a_s );

      const Tensor<1,dim> tidal_gravity = T_factor *
                                          (dTstar_gradient + dT0_gradient);

      RadialConstant<dim> radialconstant;
      return radialconstant.gravity_vector(p) + tidal_gravity;
    }


    template <int dim>
    void
    RadialWithTidalPotential<dim>::declare_parameters (ParameterHandler &prm)
    {
      RadialLinear<dim>::declare_parameters(prm);
      prm.enter_subsection("Gravity model");
      {
        prm.enter_subsection("Radial with tidal potential");
        {
          prm.declare_entry ("Mass of perturbing body", "1.898e27",
                             Patterns::Double (),
                             "Mass of body that perturbs gravity of modeled body. "
                             "The default value is chosen for modeling Europa, therefore, it is the mass of Jupiter. "
                             "Units is $kg$.");
          prm.declare_entry ("Semimajor axis of orbit", "670900000",
                             Patterns::Double (),
                             "The length of the semimajor axis of the orbit that cause the tidal perturbation. "
                             "For example, tidal perturbation on Europa happens by Europa orbiting Jupiter, "
                             "and that on Earth, if Moon is in consideration, happens by Moon orbiting Earth. "
                             "The default value is for the semimajor axis of Europa's orbit. "
                             "Units is $m$.");
          prm.declare_entry ("Angular rate of nonsynchronous rotation", "0.036",
                             Patterns::Double (),
                             "Angular rate of nonsynchronous rotation (NSR). "
                             "This works for the modeled body having decoupled rotation between interior layers. "
                             "The default value is the angular rate of Europa's icy shell. "
                             "Units is $degrees/year$ when 'Use years instead of seconds' is true, "
                             "and $degress/second$ when 'Use years instead of seconds' is false. ");
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
          const double time_scale = this->get_parameters().convert_to_years ? constants::year_in_seconds : 1.0;
          b = prm.get_double ("Angular rate of nonsynchronous rotation") * constants::degree_to_radians / time_scale;
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
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
