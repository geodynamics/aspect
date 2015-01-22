/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/initial_conditions/harmonic_perturbation.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/utilities.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace aspect
{
  namespace InitialConditions
  {
    // NOTE: this module uses the Boost spherical harmonics package which is not designed
    // for very high order (> 100) spherical harmonics computation. If you use harmonic
    // perturbations of a high order be sure to confirm the accuracy first.
    // For more information, see:
    // http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_poly/sph_harm.html
    //
    template <int dim>
    double
    HarmonicPerturbation<dim>::
    initial_temperature (const Point<dim> &position) const
    {

      // use either the user-input reference temperature as background temperature
      // (incompressible model) or the adiabatic temperature profile (compressible model)
      const double background_temperature = this->get_material_model().is_compressible() ?
                                            this->get_adiabatic_conditions().temperature(position) :
                                            reference_temperature;

      // s = fraction of the way from
      // the inner to the outer
      // boundary; 0<=s<=1
      const double s = this->get_geometry_model().depth(position) / this->get_geometry_model().maximal_depth();

      const double depth_perturbation = std::sin(vertical_wave_number*s*numbers::PI);


      double lateral_perturbation = 0.0;

      if (const GeometryModel::SphericalShell<dim> *
          spherical_geometry_model = dynamic_cast <const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()))
        {
          // In case of spherical shell calculate spherical coordinates
          const std_cxx1x::array<double,dim> scoord = aspect::Utilities::spherical_coordinates(position);

          if (dim==2)
            {
              // Use a sine as lateral perturbation that is scaled to the opening angle of the geometry.
              // This way the perturbation is alway 0 at the model boundaries.
              const double opening_angle = spherical_geometry_model->opening_angle()*numbers::PI/180.0;
              lateral_perturbation = std::sin(lateral_wave_number_1*scoord[1]*numbers::PI/opening_angle);
            }

          else if (dim==3)
            {
              // Spherical harmonics are only defined for order <= degree
              // and degree >= 0. Verify that it is indeed.
              Assert ( std::abs(lateral_wave_number_2) <= lateral_wave_number_1,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "order <= degree."));
              Assert ( lateral_wave_number_1 >= 0,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "degree >= 0."));
              // use a spherical harmonic function as lateral perturbation
              lateral_perturbation = boost::math::spherical_harmonic_r(lateral_wave_number_1,lateral_wave_number_2,scoord[2],scoord[1]);
            }
        }
      else if (const GeometryModel::Box<dim> *
               box_geometry_model = dynamic_cast <const GeometryModel::Box<dim>*> (&this->get_geometry_model()))
        {
          // In case of Box model use a sine as lateral perturbation
          // that is scaled to the extent of the geometry.
          // This way the perturbation is alway 0 at the model borders.
          const Point<dim> extent = box_geometry_model->get_extents();

          if (dim==2)
            {
              lateral_perturbation = std::sin(lateral_wave_number_1*position(0)*numbers::PI/extent(0));
            }
          else if (dim==3)
            {
              lateral_perturbation = std::sin(lateral_wave_number_1*position(0)*numbers::PI/extent(0))
                                     * std::sin(lateral_wave_number_2*position(1)*numbers::PI/extent(1));
            }
        }
      else
        AssertThrow (false,
                     ExcMessage ("Not a valid geometry model for the initial conditions model"
                                 "harmonic perturbation."));

      return background_temperature + magnitude * depth_perturbation * lateral_perturbation;
    }

    template <int dim>
    void
    HarmonicPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Harmonic perturbation");
        {
          prm.declare_entry ("Vertical wave number", "1",
                             Patterns::Integer (),
                             "Doubled radial wave number of the harmonic perturbation. "
                             " One equals half of a sine period over the model domain. "
                             " This allows for single up-/downswings. Negative numbers "
                             " reverse the sign of the perturbation.");
          prm.declare_entry ("Lateral wave number one", "3",
                             Patterns::Integer (),
                             "Doubled first lateral wave number of the harmonic perturbation. "
                             "Equals the spherical harmonic degree in 3D spherical shells. "
                             "In all other cases one equals half of a sine period over "
                             "the model domain. This allows for single up-/downswings. "
                             "Negative numbers reverse the sign of the perturbation but are "
                             "not allowed for the spherical harmonic case.");
          prm.declare_entry ("Lateral wave number two", "2",
                             Patterns::Integer (),
                             "Doubled second lateral wave number of the harmonic perturbation. "
                             "Equals the spherical harmonic order in 3D spherical shells. "
                             "In all other cases one equals half of a sine period over "
                             "the model domain. This allows for single up-/downswings. "
                             "Negative numbers reverse the sign of the perturbation.");
          prm.declare_entry ("Magnitude", "1.0",
                             Patterns::Double (0),
                             "The magnitude of the Harmonic perturbation.");
          prm.declare_entry ("Reference temperature", "1600.0",
                             Patterns::Double (0),
                             "The reference temperature that is perturbed by the"
                             "harmonic function. Only used in incompressible models.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    HarmonicPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial conditions");
      {
        prm.enter_subsection("Harmonic perturbation");
        {
          vertical_wave_number = prm.get_integer ("Vertical wave number");
          lateral_wave_number_1 = prm.get_integer ("Lateral wave number one");
          lateral_wave_number_2 = prm.get_integer ("Lateral wave number two");
          magnitude = prm.get_double ("Magnitude");
          reference_temperature = prm.get_double ("Reference temperature");
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
  namespace InitialConditions
  {
    ASPECT_REGISTER_INITIAL_CONDITIONS(HarmonicPerturbation,
                                       "harmonic perturbation",
                                       "An initial temperature field in which the temperature "
                                       "is perturbed following a harmonic function (spherical "
                                       "harmonic or sine depending on geometry and dimension) "
                                       "in lateral and radial direction from an otherwise "
                                       "constant temperature (incompressible model) or adiabatic "
                                       "reference profile (compressible model).")
  }
}
