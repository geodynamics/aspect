/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/quadrature_points.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      QuadraturePoints<dim>::generate_particles(Particles::ParticleHandler<dim> &particle_handler)
      {
        const Quadrature<dim> &quadrature_formula
          = this->introspection().quadratures.velocities;

        Particles::Generators::regular_reference_locations(
          this->get_triangulation(),
          quadrature_formula.get_points(),
          particle_handler,
          this->get_mapping());
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
      ASPECT_REGISTER_PARTICLE_GENERATOR(QuadraturePoints,
                                         "quadrature points",
                                         "Generates particles at the quadrature points of each active cell of "
                                         "the triangulation. Here, Gauss quadrature of degree (velocity\\_degree + 1), "
                                         "is used similarly to the assembly of Stokes matrix.")
    }
  }
}
