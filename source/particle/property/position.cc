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

#include <aspect/particle/property/position.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      Position<dim>::initialize_particle(std::vector<double> &data,
                                         const Point<dim> &position,
                                         const Vector<double> &)
      {
        for (unsigned int i = 0; i < dim; ++i)
          data.push_back(position[i]);
      }

      template <int dim>
      void
      Position<dim>::update_particle(unsigned int data_position,
                                   std::vector<double> &data,
                                   const Point<dim> &position,
                                   const Vector<double> &)
      {
        for (unsigned int i = 0; i < dim; ++i)
          data[data_position++] = position[i];
      }

      template <int dim>
      bool
      Position<dim>::need_update()
      {
        return true;
      }

      template <int dim>
      unsigned int
      Position<dim>::data_len() const
      {
        return dim;
      }

      /**
       * Set up the MPI data type information for the Position type
       *
       * @param [in,out] data_info Vector to append MPIDataInfo objects to
       */
      template <int dim>
      void
      Position<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info) const
      {
        for (unsigned int i = 0; i < data_len(); i++)
          {
            std::ostringstream field_name;
            field_name << "position_" << i;
            data_info.push_back(aspect::Particle::MPIDataInfo(field_name.str(), 1));
          }
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(Position,
                                        "position",
                                        "Implementation of a plugin in which the tracer "
                                        "property is defined as the recent position. "
                                        "This is used for example for "
                                        "multi-step particle integrators, because at "
                                        "the time of particle movement this property is "
                                        "not updated yet, so it can be used as the old "
                                        "position in the integrator.\n\n")
    }
  }
}

