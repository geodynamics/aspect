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

#include <aspect/particle/generator/uniform_radial.h>

#include <aspect/utilities.h>
#include <deal.II/grid/grid_tools.h>

#include <boost/random.hpp>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      // Generate a uniform radial distribution of particles over the specified region
      // in the computational domain

          /**
           * Constructor.
           *
           * @param[in] The MPI communicator for synchronizing particle generation.
           */
        template <int dim>
        UniformRadial<dim>::UniformRadial() {}

          /**
           * Generate a uniformly randomly distributed set of particles in the current triangulation.
           */
          // TODO: fix the particle system so it works even with processors assigned 0 cells
        template <int dim>
        void
		UniformRadial<dim>::generate_particles(Particle::World<dim> &world,
                                                        const double total_num_particles)
        {
           unsigned int startID = 0;

           // Create the array of shell to deal with
           std::vector<double> shell_radius(radial_layers);
           const double shell_seperation = (P_max[0] - P_min[0]) / (radial_layers-1);
           std::vector<unsigned int> ppr(radial_layers);
           double radiusTotal = 0;

           for (unsigned int i = 0; i < radial_layers; ++i)
             {
               // Calculate radius of each shell
               radiusTotal += shell_radius[i] = P_min[0] + (shell_seperation * i);
             }

           for (unsigned int i = 0; i < radial_layers; ++i)
             {
               // Calculate amount of particles per shell.
               // Number of particles depend on the portion of the radius that this shell is in (i.e., more radius = more particles)
               ppr[i] = round(total_num_particles * shell_radius[i] / radiusTotal);
             }

           uniform_radial_particles_in_subdomain(world, startID, shell_radius, ppr);

        }

        template <int dim>
        void
        UniformRadial<dim>::uniform_radial_particles_in_subdomain (Particle::World<dim> &world,
                                                    const unsigned int start_id,
                                                    const std::vector<double> &shell_radius,
                                                    const std::vector<unsigned int> &particlesPerRadius)
        {
          unsigned int cur_id;
          cur_id = start_id;
          std_cxx11::array<double,dim> spherical_coordinates;

          if (dim == 3)
            {
              for (unsigned int i = 0; i < radial_layers; i++)
                {
                  spherical_coordinates[0] = shell_radius[i];
                  const int thetaParticles = floor(sqrt(particlesPerRadius[i]));
                  const int phiParticles = particlesPerRadius[i] / thetaParticles;
                  const double phiSeperation = (P_max[1] - P_min[1]) / phiParticles;

                  std::vector<int> ppPh(phiParticles);
                  int j = 0;

                  for (double phi = P_min[1]; phi < P_max[1]; phi += phiSeperation, j++)
                    {
                      //Average value of sin(n) from 0 to 180 degrees is (2/pi)
                      ppPh[j] = (thetaParticles * sin(phi / 180 * M_PI) * (M_PI / 2)) + 1;
                    }

                  j = 0;
                  for (double phi = P_min[1]; phi < P_max[1]; phi += phiSeperation, j++)
                    {
                      spherical_coordinates[1] = phi;
                      const double thetaSeperation = (P_max[2] - P_max[1]) / ppPh[j];
                      for (double theta = P_min[2]; theta < P_max[2]; theta += thetaSeperation)
                        {
                          spherical_coordinates[2] = theta;
                          const Point<dim> newPoint = Utilities::cartesian_coordinates<dim>(spherical_coordinates) + P_center;

                          cur_id++;

                          typename parallel::distributed::Triangulation<dim>::active_cell_iterator it =
                              (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), newPoint)).first;

                          if (it->is_locally_owned())
                            {
                              //Only try to add the point if the cell it is in, is on this processor
                              BaseParticle<dim> new_particle(newPoint, cur_id);
                              world.add_particle(new_particle, std::make_pair(it->level(), it->index()));
                            }
                        }
                    }
                }
            }
          else if (dim == 2)
            {
              for (unsigned int i = 0; i < radial_layers; i++)
                {
                  spherical_coordinates[0] = shell_radius[i];
                  const double phiSeperation = (P_max[1] - P_min[1]) / particlesPerRadius[i];

                  for (double phi = P_min[1]; phi < P_max[1]; phi += phiSeperation)
                    {
                      spherical_coordinates[1] = phi;
                      const Point<dim> newPoint = Utilities::cartesian_coordinates<dim>(spherical_coordinates) + P_center;

                      //Modify the find_active_cell_around_point to only search for nearest vertex  and adj. cells, instead of searching all cells in the simulation
                      typename parallel::distributed::Triangulation<dim>::active_cell_iterator it =
                          (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), newPoint)).first;

                      if (it->is_locally_owned())
                        {
                          //Only try to add the point if the cell it is in, is on this processor
                          BaseParticle<dim> new_particle(newPoint, cur_id++);
                          world.add_particle(new_particle, std::make_pair(it->level(), it->index()));
                        }
                    }
                }
            }

        }


        template <int dim>
        void
        UniformRadial<dim>::declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Postprocess");
          {
            prm.enter_subsection("Tracers");
            {
              prm.enter_subsection("Generator");
              {
                prm.enter_subsection("Uniform radial");
                {
                prm.declare_entry ("Minimal radius", "0",
                                   Patterns::Double (),
                                   "Minimal radial coordinate for the region of tracers.");
                prm.declare_entry ("Maximal radius", "1",
                                   Patterns::Double (),
                                   "Maximal radial coordinate for the region of tracers.");
                prm.declare_entry ("Minimal longitude", "0",
                                   Patterns::Double (),
                                   "Minimal longitude coordinate for the region of tracers.");
                prm.declare_entry ("Maximal longitude", "3.1415",
                                   Patterns::Double (),
                                   "Maximal longitude coordinate for the region of tracers.");
                prm.declare_entry ("Minimal latitude", "0",
                                   Patterns::Double (),
                                   "Minimal latitude coordinate for the region of tracers.");
                prm.declare_entry ("Maximal latitude", "3.1415",
                                   Patterns::Double (),
                                   "Maximal latitude coordinate for the region of tracers.");
                prm.declare_entry ("Radial layers", "1",
                                   Patterns::Integer(1),
                                   "The number of radial layers of particles that will be generated.");
                prm.declare_entry ("Center x", "0",
                                   Patterns::Double (),
                                   "x coordinate for the center of the region of tracers.");
                prm.declare_entry ("Center y", "0",
                                   Patterns::Double (),
                                   "y coordinate for the center of the region of tracers.");
                prm.declare_entry ("Center z", "0",
                                   Patterns::Double (),
                                   "z coordinate for the center of the region of tracers.");
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }


        template <int dim>
        void
        UniformRadial<dim>::parse_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Postprocess");
          {
            prm.enter_subsection("Tracers");
            {
              prm.enter_subsection("Generator");
              {
                prm.enter_subsection("Uniform radial");
                {
                  P_min[0] = prm.get_double ("Minimal radius");
                  P_max[0] = prm.get_double ("Maximal radius");
                  P_min[1] = prm.get_double ("Minimal longitude");
                  P_max[1] = prm.get_double ("Maximal longitude");

                  P_center[0] = prm.get_double ("Center x");
                  P_center[1] = prm.get_double ("Center y");

                  if (dim ==3)
                    {
                      P_min[2] = prm.get_double ("Minimal latitude");
                      P_max[2] = prm.get_double ("Maximal latitude");
                      P_center[2] = prm.get_double ("Center z");
                    }

                  radial_layers = prm.get_integer("Radial layers");
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
    ASPECT_REGISTER_PARTICLE_GENERATOR(UniformRadial,
                                               "uniform radial",
                                               "Generate a uniform distribution of particles"
                                               "over a spherical domain in 2D or 3D. Note that in order "
                                               "to produce a regular distribution the number of generated "
                                               "tracers might not exactly match the one specified in the "
                                               "input file.")
    }
  }
}
