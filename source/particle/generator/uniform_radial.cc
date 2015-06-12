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


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      UniformRadial<dim>::UniformRadial() {}

      template <int dim>
      void
      UniformRadial<dim>::generate_particles(World<dim> &world)
      {
        // Create the array of shell to deal with
        const double radius_spacing = (P_max[0] - P_min[0]) / fmax(radial_layers-1,1);

        // Calculate amount of particles per shell.
        // The number of particles depend on the fraction of the area
        // (or length in 2D) that this shell occupies compared to the total domain
        std::vector<unsigned int> particles_per_radius(radial_layers);
        if (dim == 2)
          {
            double total_radius = 0;
            for (unsigned int i = 0; i < radial_layers; ++i)
              total_radius += P_min[0] + (radius_spacing * i);
            for (unsigned int i = 0; i < radial_layers; ++i)
              {
                const double radius = P_min[0] + (radius_spacing * i);
                particles_per_radius[i] = round(n_tracers * radius / total_radius);
              }
          }
        else if (dim == 3)
          {
            double total_area = 0;
            for (unsigned int i = 0; i < radial_layers; ++i)
              total_area += std::pow(P_min[0] + (radius_spacing * i),2);
            for (unsigned int i = 0; i < radial_layers; ++i)
              {
                const double area = std::pow(P_min[0] + (radius_spacing * i),2);
                particles_per_radius[i] = round(n_tracers * area / total_area);
              }
          }
        else
          ExcNotImplemented();

        unsigned int cur_id = 0;
        std_cxx11::array<double,dim> spherical_coordinates;

        for (unsigned int i = 0; i < radial_layers; ++i)
          {
            spherical_coordinates[0] = P_min[0] + (radius_spacing * i);
            if (dim == 2)
              {
                const double phi_spacing = (P_max[1] - P_min[1]) / fmax(particles_per_radius[i]-1,1);

                for (unsigned int j = 0; j < particles_per_radius[i]; ++j)
                  {
                    spherical_coordinates[1] = P_min[1] + j * phi_spacing;
                    const Point<dim> newPoint = Utilities::cartesian_coordinates<dim>(spherical_coordinates) + P_center;
                    generate_particle(newPoint,cur_id++,world);
                  }
              }
            else if (dim == 3)
              {
                const unsigned int theta_particles = round(sqrt(particles_per_radius[i]));
                const unsigned int phi_particles = round ((double) particles_per_radius[i] / (double) theta_particles);
                const double theta_spacing = (P_max[2] - P_min[2]) / fmax(theta_particles-1,1);

                for (unsigned int j = 0; j < theta_particles; ++j)
                  {
                    spherical_coordinates[2] = P_min[2] + j * theta_spacing;

                    //Average value of sin(n) from 0 to 180 degrees is (2/pi)
                    const int adjusted_phi_particles = std::max(static_cast<int> (phi_particles * std::sin(spherical_coordinates[2])),1);
                    const double phi_spacing = (P_max[1] - P_min[1]) / fmax(adjusted_phi_particles-1,1);
                    for (unsigned int k = 0; k < adjusted_phi_particles; ++k)
                      {
                        spherical_coordinates[1] = P_min[1] + k * phi_spacing;
                        const Point<dim> newPoint = Utilities::cartesian_coordinates<dim>(spherical_coordinates) + P_center;
                        generate_particle(newPoint,cur_id++,world);
                      }
                  }
              }
            else
              ExcNotImplemented();
          }
      }


      template <int dim>
      void
      UniformRadial<dim>::generate_particle(const Point<dim> &position,
                                            const unsigned int id,
                                            World<dim> &world)
      {
        const typename parallel::distributed::Triangulation<dim>::active_cell_iterator it =
          (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), position)).first;;

        if (it->is_locally_owned())
          {
            //Only try to add the point if the cell it is in, is on this processor
            Particle<dim> new_particle(position, id);
            world.add_particle(new_particle, std::make_pair(it->level(), it->index()));
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
            prm.declare_entry ("Number of tracers", "1000",
                               Patterns::Double (0),
                               "Total number of tracers to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

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
            n_tracers    = static_cast<unsigned int>(prm.get_double ("Number of tracers"));

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
