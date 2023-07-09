/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/random_gaussian_perturbation.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/simulator_access.h>

#include <random>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    void
    RandomGaussianPerturbation<dim>::initialize ()
    {
      // Get dimension of the model from the geometry model.
      Point<dim> min_coordinates;
      Point<dim> max_coordinates;

      perturbation_centers.resize(n_perturbations);
      perturbation_magnitudes.resize(n_perturbations);

      if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model()))
        {
          const GeometryModel::SphericalShell<dim> &spherical_geometry_model =
            Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());

          for (unsigned int d=0; d<dim; ++d)
            {
              min_coordinates[d] = -spherical_geometry_model.outer_radius();
              max_coordinates[d] = spherical_geometry_model.outer_radius();
            }
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>> (this->get_geometry_model()))
        {
          const GeometryModel::Sphere<dim> &sphere_geometry_model =
            Plugins::get_plugin_as_type<const GeometryModel::Sphere<dim>> (this->get_geometry_model());

          for (unsigned int d=0; d<dim; ++d)
            {
              min_coordinates[d] = -sphere_geometry_model.radius();
              max_coordinates[d] = sphere_geometry_model.radius();
            }
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>> (this->get_geometry_model()))
        {
          const GeometryModel::Chunk<dim> &chunk_geometry_model =
            Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>> (this->get_geometry_model());

          for (unsigned int d=0; d<dim; ++d)
            {
              min_coordinates[d] = -chunk_geometry_model.outer_radius();
              max_coordinates[d] = chunk_geometry_model.outer_radius();
            }
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model()))
        {
          const GeometryModel::Box<dim> &box_geometry_model =
            Plugins::get_plugin_as_type<const GeometryModel::Box<dim>> (this->get_geometry_model());

          min_coordinates = box_geometry_model.get_origin();
          max_coordinates = min_coordinates + box_geometry_model.get_extents();
        }
      else if (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
        {
          const GeometryModel::TwoMergedBoxes<dim> &two_merged_boxes_geometry_model =
            Plugins::get_plugin_as_type<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model());

          min_coordinates = two_merged_boxes_geometry_model.get_origin();
          max_coordinates = min_coordinates + two_merged_boxes_geometry_model.get_extents();
        }
      else
        AssertThrow (false,
                     ExcMessage ("Not a valid geometry model for the initial conditions model "
                                 "'random Gaussian perturbation'."));

      // Initialize random locations of perturbations and check they are in the model domain.
      // Use a fixed number as seed for random generator.
      // This is important if we run the code on more than 1 processor.
      std::mt19937 generator(1);
      std::uniform_real_distribution<double> random_location(0.0,1.0);
      std::uniform_real_distribution<double> random_magnitude(-max_magnitude, max_magnitude);

      for (unsigned int n=0; n<n_perturbations; ++n)
        {
          do
            {
              for (unsigned int d=0; d<dim; ++d)
                perturbation_centers[n][d] = min_coordinates[d] + random_location(generator) * (max_coordinates[d] - min_coordinates[d]);
            }
          while (!this->get_geometry_model().point_is_in_domain(perturbation_centers[n]));

          perturbation_magnitudes[n] = random_magnitude(generator);
        }
    }



    template <int dim>
    double
    RandomGaussianPerturbation<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      double temperature_perturbation = 0;
      for (unsigned int n=0; n<n_perturbations; ++n)
        {
          const double distance_square = position.distance_square(perturbation_centers[n]);
          const double exponent = distance_square / (2.*width*width);

          temperature_perturbation += perturbation_magnitudes[n] * std::exp(-exponent);
        }

      return temperature_perturbation;
    }



    template <int dim>
    void
    RandomGaussianPerturbation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Random Gaussian perturbation");
        {
          prm.declare_entry ("Number of perturbations", "100",
                             Patterns::Integer (),
                             "Total number of perturbations to be introduced into the model. "
                             "Perturbations will be placed at random locations within the "
                             "model domain.");
          prm.declare_entry ("Maximum magnitude", "25.0",
                             Patterns::Double (0.),
                             "The maximum magnitude of the Gaussian perturbation. For each "
                             "perturbation, a random magnitude between plus and minus the "
                             "maximum magnitude will be chosen. "
                             "Units: \\si{\\kelvin}.");
          prm.declare_entry ("Width", "1000.0",
                             Patterns::Double (0.),
                             "The Gaussian RMS width of the perturbations. "
                             "Units: \\si{\\meter}.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    RandomGaussianPerturbation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Random Gaussian perturbation");
        {
          n_perturbations = prm.get_integer ("Number of perturbations");
          max_magnitude = prm.get_double ("Maximum magnitude");
          width = prm.get_double ("Width");
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
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(RandomGaussianPerturbation,
                                              "random Gaussian perturbation",
                                              "An initial temperature field in which the temperature "
                                              "is perturbed from a temperature of zero following a "
                                              "given number of Gaussian perturbations placed randomly "
                                              "throughout the model domain. The number, width, and "
                                              "maximum magnitude of the perturbations can be chosen "
                                              "as model parameters. "
                                              "This plugin is meant to be used in combination with "
                                              "another initial temperature model that determines the "
                                              "background temperature (such as the 'function' or the "
                                              "'adiabatic' plugin) using the 'add' operator to combine "
                                              "them.")
  }
}
