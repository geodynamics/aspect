/*
  Copyright (C) 2018 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPEC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/postprocess/gravity_point_values.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/global.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    GravityPointValues<dim>::execute (TableHandler &)
    {
      // Get quadrature formula and increase the degree of quadrature over the velocity
      // element degree.
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+quadrature_degree_increase);
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      // Get the value of the outer radius and inner radius
      double model_outer_radius;
      double model_inner_radius;
      if (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != nullptr)
        {
          model_outer_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                               (this->get_geometry_model()).outer_radius();
          model_inner_radius = dynamic_cast<const GeometryModel::SphericalShell<dim>&>
                               (this->get_geometry_model()).inner_radius();
        }
      else if (dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != nullptr)
        {
          model_outer_radius = dynamic_cast<const GeometryModel::Chunk<dim>&>
                               (this->get_geometry_model()).outer_radius();
          model_inner_radius = dynamic_cast<const GeometryModel::Chunk<dim>&>
                               (this->get_geometry_model()).inner_radius();
        }
      else if (dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != nullptr)
        {
          model_outer_radius = dynamic_cast<const GeometryModel::Sphere<dim>&>
                               (this->get_geometry_model()).radius();
          model_inner_radius = 0;
        }
      else
        {
          Assert (false, ExcMessage ("This initial condition can only be used if the geometry "
                                     "is a sphere, a spherical shell, a chunk or an "
                                     "ellipsoidal chunk."));
          model_outer_radius = 1;
          model_inner_radius = 0;
        }

      // Get the value of the universal gravitational constant
      const double G = aspect::constants::big_g;

      // now write all data to the file of choice. start with a pre-amble that
      // explains the meaning of the various fields
      const std::string filename = (this->get_output_directory() +
                                    "output_gravity");
      std::ofstream output (filename.c_str());
      output << "# 1: position_satellite_r" << '\n'
             << "# 2: position_satellite_phi" << '\n'
             << "# 3: position_satellite_theta" << '\n'
             << "# 4: position_satellite_x" << '\n'
             << "# 5: position_satellite_y" << '\n'
             << "# 6: position_satellite_z" << '\n'
             << "# 7: gravity_x" << '\n'
             << "# 8: gravity_y" << '\n'
             << "# 9: gravity_z" << '\n'
             << "# 10: gravity_norm" << '\n'
             << "# 11: gravity_theory" << '\n'
             << "# 12: gravity potential" << '\n'
             << "# 13: gravity_anomaly_x" << '\n'
             << "# 14: gravity_anomaly_y" << '\n'
             << "# 15: gravity_anomaly_z" << '\n'
             << "# 16: gravity_anomaly_norm" << '\n'
             << '\n';

      // Storing cartesian coordinate, density and JxW at local quadrature points in a vector
      // avoids to use MaterialModel and fe_values within the loops. Because postprocessor
      // run in parallel, the total number of local quadrature points has to be determined:
      const unsigned int n_locally_owned_cells = (this->get_triangulation().n_locally_owned_active_cells());
      const unsigned int n_quadrature_points_per_cell = quadrature_formula.size();

      // Declare the vector 'density_JxW' to store the density at quadrature points. The
      // density and the JxW are here together for simplicity in the equation (both variables
      // only appear together):
      std::vector<double> density_JxW (n_locally_owned_cells * n_quadrature_points_per_cell);
      std::vector<double> relative_density_JxW (n_locally_owned_cells * n_quadrature_points_per_cell);

      // Declare the vector 'position_point' to store the position of quadrature points:
      std::vector<Point<dim> > position_point (n_locally_owned_cells * n_quadrature_points_per_cell);

      // The following loop perform the storage of the position and density * JxW values
      // at local quadrature points:
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      MaterialModel::MaterialModelInputs<dim> in(quadrature_formula.size(),this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature_formula.size(),this->n_compositional_fields());
      unsigned int local_cell_number = 0;
      for (; cell!=endc; ++cell)
        {
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              const std::vector<Point<dim> > &position_point_cell = fe_values.get_quadrature_points();
              in.reinit(fe_values, cell, this->introspection(), this->get_solution(), false);
              this->get_material_model().evaluate(in, out);
              for (unsigned int q = 0; q < n_quadrature_points_per_cell; ++q)
                {
                  density_JxW[local_cell_number * n_quadrature_points_per_cell + q] = out.densities[q] * fe_values.JxW(q);
                  relative_density_JxW[local_cell_number * n_quadrature_points_per_cell + q] = (out.densities[q]-reference_density) * fe_values.JxW(q);
                  position_point[local_cell_number * n_quadrature_points_per_cell + q] = position_point_cell[q];
                }
              ++local_cell_number;
            }
        }

      // This is the main loop which computes gravity acceleration and potential at a
      // point located at the spherical coordinate [r, phi, theta]:
      // loop on the radius - satellite position [r, , ]
      for (unsigned int h=0; h < number_points_radius; ++h)
        {
          std::array<double,dim> satellite_coordinate;
          satellite_coordinate[0] = minimum_radius + ((maximum_radius - minimum_radius) / number_points_radius) * h;

          // loop on the longitude - satellite position [ , phi, ] in radian:
          for (unsigned int i=0; i < number_points_longitude; ++i)
            {
              satellite_coordinate[1] = (minimum_longitude + ((maximum_longitude - minimum_longitude) / number_points_longitude) * i) * numbers::PI / 180.;

              // loop on latitude - satllite position [ , , theta] in radian:
              for (unsigned int j=0; j < number_points_latitude; ++j)
                {
                  satellite_coordinate[2] = (minimum_latitude + ((maximum_latitude - minimum_latitude) / number_points_latitude) * j) * numbers::PI / 180.;

                  // The spherical coordinates are shifted into cartesian to allow simplification in the mathematical equation.
                  const Point<dim> position_satellite = Utilities::Coordinates::spherical_to_cartesian_coordinates<dim>(satellite_coordinate);

                  // For each point (i.e. satellite), the fourth integral goes over cells and quadrature
                  // points to get the unique distance between those, indispensable to calculate
                  // gravity vector components x,y,z (in tensor), and potential.
                  Tensor<1,dim> local_g;
                  Tensor<1,dim> local_g_anomaly;
                  double local_g_potential = 0;
                  cell = this->get_dof_handler().begin_active();
                  local_cell_number = 0;
                  for (; cell!=endc; ++cell)
                    {
                      if (cell->is_locally_owned())
                        {
                          for (unsigned int q = 0; q < n_quadrature_points_per_cell; ++q)
                            {
                              const double dist = (position_satellite - position_point[local_cell_number * n_quadrature_points_per_cell + q]).norm();
                              const double KK = G * density_JxW[local_cell_number * n_quadrature_points_per_cell + q] / pow(dist,3);
                              const double relative_KK = G * relative_density_JxW[local_cell_number * n_quadrature_points_per_cell + q] / pow(dist,3);
                              local_g += KK * (position_satellite - position_point[local_cell_number * n_quadrature_points_per_cell + q]);
                              local_g_anomaly += relative_KK * (position_satellite - position_point[local_cell_number * n_quadrature_points_per_cell + q]);
                              local_g_potential -= G * density_JxW[local_cell_number * n_quadrature_points_per_cell + q] / dist;
                            }
                          ++local_cell_number;
                        }
                    }

                  // Sum local gravity components over global domain
                  const Tensor<1,dim> g
                    = Utilities::MPI::sum (local_g, this->get_mpi_communicator());
                  const double g_potential
                    = Utilities::MPI::sum (local_g_potential, this->get_mpi_communicator());
                  const Tensor<1,dim> g_anomaly
                    = Utilities::MPI::sum (local_g_anomaly, this->get_mpi_communicator());

                  // analytical solution to calculate the theoretical gravity from a uniform density model.
                  // can only be used if concentric density profile
                  double g_theory = 0;
                  if (satellite_coordinate[0] <= model_inner_radius)
                    {
                      g_theory = 0;
                    }
                  else if ((satellite_coordinate[0] > model_inner_radius) && (satellite_coordinate[0] < model_outer_radius))
                    {
                      g_theory = (G * numbers::PI * 4/3 * reference_density * (satellite_coordinate[0] -
                                                                               (pow(model_inner_radius,3) / pow(satellite_coordinate[0],2))));
                    }
                  else
                    {
                      g_theory = (G * numbers::PI * 4/3 * reference_density * (pow(model_outer_radius,3) -
                                                                               pow(model_inner_radius,3)) / pow(satellite_coordinate[0],2));
                    }

                  // write output
                  if (dealii::Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
                    {
                      output << satellite_coordinate[0] << ' '
                             << satellite_coordinate[1] *180. / numbers::PI << ' '
                             << satellite_coordinate[2] *180. / numbers::PI << ' '
                             << position_satellite[0] << ' '
                             << position_satellite[1] << ' '
                             << position_satellite[2] << ' '
                             << g << ' '
                             << g.norm() << ' '
                             << g_theory << ' '
                             << g_potential
                             << g_anomaly << ' '
                             << g_anomaly.norm() << ' '
                             << '\n';
                    }
                }
            }
        }
      return std::pair<std::string, std::string> ("gravity computation file:",
                                                  filename);
    }

    template <int dim>
    void
    GravityPointValues<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          prm.declare_entry ("Quadrature degree increase", "1",
                             Patterns::Double (0.0),
                             "Quadrature degree increase over the velocity element "
                             "degree may be required when gravity is calculated near "
                             "the surface or inside the model. An increase in the "
                             "quadrature element adds accuracy to the gravity "
                             "solution from noise due to the model grid.");
          prm.declare_entry ("Number points radius", "1",
                             Patterns::Double (0.0),
                             "Gravity may be calculated for a set of points along "
                             "the radius (e.g. depth profile) between a minimum and "
                             "maximum radius.");
          prm.declare_entry ("Number points longitude", "1",
                             Patterns::Double (0.0),
                             "Gravity may be calculated for a sets of points along "
                             "the longitude (e.g. gravity map) between a minimum and "
                             "maximum longitude.");
          prm.declare_entry ("Number points latitude", "1",
                             Patterns::Double (0.0),
                             "Gravity may be calculated for a sets of points along "
                             "the latitude (e.g. gravity map) between a minimum and "
                             "maximum latitude.");
          prm.declare_entry ("Minimum radius", "0",
                             Patterns::Double (0.0),
                             "Minimum radius may be defined in or outside the model.");
          prm.declare_entry ("Maximum radius", "0",
                             Patterns::Double (0.0),
                             "Maximum radius can be defined in or outside the model.");
          prm.declare_entry ("Minimum longitude", "0",
                             Patterns::Double (0.0,360.0),
                             "Gravity may be calculated for a sets of points along "
                             "the longitude between a minimum and maximum longitude.");
          prm.declare_entry ("Minimum latitude", "0",
                             Patterns::Double (0.0,180.0),
                             "Gravity may be calculated for a sets of points along "
                             "the latitude between a minimum and maximum latitude.");
          prm.declare_entry ("Maximum longitude", "360",
                             Patterns::Double (0.0,360.0),
                             "Gravity may be calculated for a sets of points along "
                             "the longitude between a minimum and maximum longitude.");
          prm.declare_entry ("Maximum latitude", "180",
                             Patterns::Double (0.0,180.0),
                             "Gravity may be calculated for a sets of points along "
                             "the latitude between a minimum and maximum latitude.");
          prm.declare_entry ("Reference density", "3300",
                             Patterns::Double (0.0),
                             "Gravity anomalies are computed using densities anomalies "
                             "relative to a reference density.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    GravityPointValues<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("Gravity calculation");
        {
          quadrature_degree_increase = prm.get_double ("Quadrature degree increase");
          number_points_radius    = prm.get_double ("Number points radius");
          number_points_longitude = prm.get_double ("Number points longitude");
          number_points_latitude  = prm.get_double ("Number points latitude");
          minimum_radius    = prm.get_double ("Minimum radius");
          maximum_radius    = prm.get_double ("Maximum radius");
          minimum_longitude = prm.get_double ("Minimum longitude");
          maximum_longitude = prm.get_double ("Maximum longitude");
          minimum_latitude  = prm.get_double ("Minimum latitude");
          maximum_latitude  = prm.get_double ("Maximum latitude");
          reference_density  = prm.get_double ("Reference density");
          AssertThrow (dynamic_cast<const GeometryModel::SphericalShell<dim>*> (&this->get_geometry_model()) != nullptr ||
                       dynamic_cast<const GeometryModel::Sphere<dim>*> (&this->get_geometry_model()) != nullptr ||
                       dynamic_cast<const GeometryModel::Chunk<dim>*> (&this->get_geometry_model()) != nullptr,
                       ExcMessage ("This postprocessor can only be used if the geometry "
                                   "is a sphere, spherical shell or spherical chunk."));
          AssertThrow (! this->get_material_model().is_compressible(),
                       ExcMessage("This postprocessor was only tested for incompressible models so far."));
          AssertThrow (dim==3,
                       ExcMessage("This postprocessor was only tested for 3D models."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(GravityPointValues,
                                  "gravity calculation",
                                  "A postprocessor that computes gravity and gravity potential "
                                  "for a set of points (e.g. satellites) in or above the model "
                                  "surface for a user-defined range of latitudes, longitudes "
                                  "and radius.")
  }
}
