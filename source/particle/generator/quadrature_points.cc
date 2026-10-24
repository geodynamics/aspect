// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
