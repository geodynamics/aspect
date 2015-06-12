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

#include <aspect/particle/generator/uniform_box.h>

#include <boost/random.hpp>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/grid/grid_tools.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      UniformBox<dim>::UniformBox() {}

      template <int dim>
      void
      UniformBox<dim>::generate_particles(World<dim> &world)
      {
        unsigned int cur_id = 0;
        const Tensor<1,dim> P_diff = P_max - P_min;

        double volume(1.0);
        for (unsigned int i = 0; i < dim; ++i)
          volume *= P_diff[i];

        std_cxx11::array<double,dim> nParticles;
        std_cxx11::array<double,dim> spacing;

        // Calculate separation of particles
        for (unsigned int i = 0; i < dim; ++i)
          {
            nParticles[i] = round(std::pow(n_tracers * std::pow(P_diff[i],dim) / volume, 1./dim));
            spacing[i] = P_diff[i] / fmax(nParticles[i] - 1,1);
          }

        for (unsigned int i = 0; i < nParticles[0]; ++i)
          {
            for (unsigned int j = 0; j < nParticles[1]; ++j)
              {
                if (dim == 2)
                  generate_particle(Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1]),cur_id++,world);
                else if (dim == 3)
                  for (unsigned int k = 0; k < nParticles[2]; ++k)
                    generate_particle(Point<dim> (P_min[0]+i*spacing[0],P_min[1]+j*spacing[1],P_min[2]+k*spacing[2]),cur_id++,world);
                else
                  ExcNotImplemented();
              }
          }
      }


      template <int dim>
      void
      UniformBox<dim>::generate_particle(const Point<dim> &position,
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
      UniformBox<dim>::declare_parameters (ParameterHandler &prm)
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
              prm.enter_subsection("Uniform box");
              {
                prm.declare_entry ("Minimal x", "0",
                                   Patterns::Double (),
                                   "Minimal x coordinate for the region of tracers.");
                prm.declare_entry ("Maximal x", "1",
                                   Patterns::Double (),
                                   "Maximal x coordinate for the region of tracers.");
                prm.declare_entry ("Minimal y", "0",
                                   Patterns::Double (),
                                   "Minimal y coordinate for the region of tracers.");
                prm.declare_entry ("Maximal y", "1",
                                   Patterns::Double (),
                                   "Maximal y coordinate for the region of tracers.");
                prm.declare_entry ("Minimal z", "0",
                                   Patterns::Double (),
                                   "Minimal z coordinate for the region of tracers.");
                prm.declare_entry ("Maximal z", "1",
                                   Patterns::Double (),
                                   "Maximal z coordinate for the region of tracers.");
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
      UniformBox<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Tracers");
          {
            n_tracers    = static_cast<unsigned int>(prm.get_double ("Number of tracers"));

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Uniform box");
              {
                P_min(0) = prm.get_double ("Minimal x");
                P_max(0) = prm.get_double ("Maximal x");
                P_min(1) = prm.get_double ("Minimal y");
                P_max(1) = prm.get_double ("Maximal y");

                if (dim == 3)
                  {
                    P_min(2) = prm.get_double ("Minimal z");
                    P_max(2) = prm.get_double ("Maximal z");
                  }
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(UniformBox,
                                         "uniform box",
                                         "Generate a uniform distribution of particles "
                                         "over a rectangular domain in 2D or 3D. Note that in order "
                                         "to produce a regular distribution the number of generated "
                                         "tracers might not exactly match the one specified in the "
                                         "input file.")
    }
  }
}
