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

#include <aspect/global.h>
#include <aspect/particle/property/initial_position.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      InitialPosition<dim>::initialize_particle(std::vector<double> &data,
                                                const Point<dim> &position,
                                                const Vector<double> &,
                                                const std::vector<Tensor<1,dim> > &)
      {
        for (unsigned int i = 0; i < dim; i++)
          data.push_back(position[i]);
      }

      template <int dim>
      void
      InitialPosition<dim>::data_length(std::vector<unsigned int> &length) const
      {
        length.push_back(dim);
      }

      template <int dim>
      void
      InitialPosition<dim>::data_names(std::vector<std::string> &names) const
      {
        names.push_back("initial position");
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(InitialPosition,
                                        "initial position",
                                        "Implementation of a plugin in which the tracer "
                                        "property is given as the initial position "
                                        "of the tracer."
                                        "\n\n")
    }
  }
}

