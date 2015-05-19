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

#include <aspect/particle/property/pT_path.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      PTPath<dim>::initialize_particle(std::vector<double> &data,
                                         const Point<dim> &position,
                                         const Vector<double> &solution)
      {
          data.push_back(solution[this->introspection().component_indices.pressure]);
          data.push_back(solution[this->introspection().component_indices.temperature]);
      }

      template <int dim>
      void
      PTPath<dim>::update_particle(unsigned int data_position,
                                   std::vector<double> &data,
                                   const Point<dim> &position,
                                   const Vector<double> &solution)
      {
          data[data_position++] = solution[this->introspection().component_indices.pressure];
          data[data_position++] = solution[this->introspection().component_indices.temperature];
      }

      template <int dim>
      bool
      PTPath<dim>::need_update()
      {
        return true;
      }

      template <int dim>
      unsigned int
      PTPath<dim>::data_len() const
      {
        return 2;
      }

      /**
       * Set up the MPI data type information for the PTPath type
       *
       * @param [in,out] data_info Vector to append MPIDataInfo objects to
       */
      template <int dim>
      void
      PTPath<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info) const
      {
        data_info.push_back(aspect::Particle::MPIDataInfo("p", 1));
        data_info.push_back(aspect::Particle::MPIDataInfo("T", 1));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(PTPath,
                                        "pT path",
                                        "Implementation of a plugin in which the tracer "
                                        "property is defined as the recent pressure and "
                                        "temperature at this position. This can be used "
                                        "to calculate pressure-temperature paths of "
                                        "material points."
                                        "\n\n")
    }
  }
}

