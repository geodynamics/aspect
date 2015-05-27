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
        for (unsigned int i = 0; i < data_len(); i++)
        data.push_back(position[i]);
      }

      template <int dim>
      unsigned int
      InitialPosition<dim>::data_len() const
      {
        return dim;
      }

      /**
       * Set up the MPI data type information for the InitialPosition type
       *
       * @param [in,out] data_info Vector to append MPIDataInfo objects to
       */
      template <int dim>
      void
      InitialPosition<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info) const
      {
        data_info.push_back(aspect::Particle::MPIDataInfo("initial position", data_len()));
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

