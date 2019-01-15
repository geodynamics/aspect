/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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


#include "harmonic_perturbation_citcoms.h"
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    double
    HarmonicPerturbationCitcomS<dim>::
    initial_temperature (const Point<dim> &position) const
    {

      // s is the depth fraction of the way from the inner to the outer boundary; 0 <= depth_fraction <= 1.
      // The depth fraction is used to calculate the depth perturbation.
      const double s = this->get_geometry_model().depth(position) / this->get_geometry_model().maximal_depth();
      const double depth_perturbation = std::sin(vertical_wave_number * s * numbers::PI);

      // Get the value of the outer radius and inner radius
      double model_outer_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                                  (this->get_geometry_model()).outer_radius();
      double model_inner_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                                  (this->get_geometry_model()).inner_radius();

      // epsilon specifies the amplitude of the nonaxisymmetric perturbation. epsilon here is only valid
      // for the cubic steady-state case. For users who want to try the dodecahedral initial condition, epsilon
      // has to be set at sqrt(14/11) (Arrial et al. 2014).
      const double epsilon = 5./7. * (1 - delta);

      // In case of spherical shell, calculate spherical coordinates
      const std_cxx11::array<double,dim> scoord = aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

      // This is different than the standard harmonic perturbation in ASPECT
      double background_temperature = model_inner_radius * (scoord[0] - model_outer_radius) / scoord[0] / (model_inner_radius - model_outer_radius);

      // Use a spherical harmonic function as lateral perturbation
      std::pair<double,double> sph_harm_vals1 = Utilities::real_spherical_harmonic(lateral_wave_number_l1, lateral_wave_number_m1, scoord[2], scoord[1]);
      std::pair<double,double> sph_harm_vals2 = Utilities::real_spherical_harmonic(lateral_wave_number_l2, lateral_wave_number_m2, scoord[2], scoord[1]);

      // For historical reasons, the initial conditions module used an unnormalized real spherical harmonic
      // (see utilities.cc). The current version of the harmonic perturbation in ASPECT denormalizes the
      // return value of the real_spherical_harmonic to keep its original behavior.
      // For this benchmark, we use the original version of the harmonic perturbation, which is the same
      // than used in Arrial et al. (2014) and Zhong et al. (2008).
      // We left in comments below the denormalization which is not required here.
      double lateral_perturbation_1 = sph_harm_vals1.first; // / (lateral_wave_number_m1 == 0 ? 1.0 : std::sqrt(2.));
      double lateral_perturbation_2 = sph_harm_vals2.first; // / (lateral_wave_number_m2 == 0 ? 1.0 : std::sqrt(2.));

      return background_temperature + magnitude * depth_perturbation * (lateral_perturbation_1 + epsilon * lateral_perturbation_2);
    }

    template <int dim>
    void
    HarmonicPerturbationCitcomS<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection ("Harmonic perturbation CitcomS");
        {
          prm.declare_entry ("Vertical wave number", "1",
                             Patterns::Integer (),
                             "The radial/depth wave number of the harmonic perturbation. "
                             "One equals half of a sine period over the model domain. "
                             "This allows for single up-/downswings. Negative numbers "
                             "reverse the sign of the perturbation.");
          prm.declare_entry ("Lateral wave number l1", "3",
                             Patterns::Integer (),
                             "The degree of the spherical harmonics in a 3D spherical shell "
                             "for the lateral wave number of the axisymmetric harmonic "
                             "perturbation.");
          prm.declare_entry ("Lateral wave number m1", "2",
                             Patterns::Integer (),
                             "The order of the spherical harmonics in a 3D spherical shell "
                             "for the lateral wave number of the axisymmetric harmonic "
                             "perturbation.");
          prm.declare_entry ("Lateral wave number l2", "3",
                             Patterns::Integer (),
                             "The degree of the spherical harmonics in a 3D spherical shell "
                             "for the lateral wave number of the nonaxisymmetric harmonic "
                             "perturbation.");
          prm.declare_entry ("Lateral wave number m2", "2",
                             Patterns::Integer (),
                             "The order of the spherical harmonics in a 3D spherical shell "
                             "for the lateral wave number of the nonaxisymmetric harmonic "
                             "perturbation.");
          prm.declare_entry ("Magnitude", "1.0",
                             Patterns::Double (0),
                             "The magnitude of the Harmonic perturbation.");
          prm.declare_entry ("Delta", "1.0",
                             Patterns::Double (0),
                             "A perturbation parameter to slowly perturb the amplitude of "
                             "the nonaxisymmetric mode. Set to 0 means the perturbation is "
                             "full; set to 1 and there is no perturbation (i.e the initial "
                             "condition temperature profile is only axissymmetric).");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void
    HarmonicPerturbationCitcomS<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection ("Harmonic perturbation CitcomS");
        {
          vertical_wave_number = prm.get_integer ("Vertical wave number");
          lateral_wave_number_l1 = prm.get_integer ("Lateral wave number l1");
          lateral_wave_number_m1 = prm.get_integer ("Lateral wave number m1");
          AssertThrow (std::abs(lateral_wave_number_m1) <= lateral_wave_number_l1,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "order <= degree."));
          AssertThrow (lateral_wave_number_l1 >= 0,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "degree >= 0."));
          lateral_wave_number_l2 = prm.get_integer ("Lateral wave number l2");
          lateral_wave_number_m2 = prm.get_integer ("Lateral wave number m2");
          AssertThrow (std::abs(lateral_wave_number_m2) <= lateral_wave_number_l2,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "order <= degree."));
          AssertThrow (lateral_wave_number_l2 >= 0,
                       ExcMessage ("Spherical harmonics can only be computed for "
                                   "degree >= 0."));
          magnitude = prm.get_double ("Magnitude");
          delta = prm.get_double ("Delta");
          AssertThrow (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != 0,
                       ExcMessage ("This setup can only be used if the geometry "
                                   "is a spherical shell."));
          AssertThrow (dim==3,
                       ExcMessage ("This setup only work for 3D models."));
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
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(HarmonicPerturbationCitcomS,
                                              "harmonic perturbation CitcomS",
                                              "This setup is a different implementation than the one "
                                              "in the main code and it is based on the Arrial et al. "
                                              "(2014) setup. It generates an initial constant "
                                              "temperature field which is perturbed following a "
                                              "spherical harmonic function in lateral and radial "
                                              "direction. This setup can only be used for a hollow "
                                              "sphere.")
  }
}
