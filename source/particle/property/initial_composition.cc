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

#include <aspect/particle/property/initial_composition.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      InitialComposition<dim>::initialize_particle(std::vector<double> &data,
                                         const Point<dim> &,
                                         const Vector<double> &solution,
                                         const std::vector<Tensor<1,dim> > &)
      {
        for (unsigned int i = 0; i < data_len(); i++)
        data.push_back(solution[this->introspection().component_indices.compositional_fields[i]]);
      }

      template <int dim>
      unsigned int
      InitialComposition<dim>::data_len() const
      {
        return this->n_compositional_fields();
      }

      /**
       * Set up the MPI data type information for the InitialComposition type
       *
       * @param [in,out] data_info Vector to append MPIDataInfo objects to
       */
      template <int dim>
      void
      InitialComposition<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info) const
      {
        for (unsigned int i = 0; i < data_len(); i++)
          {
            std::ostringstream field_name;
            field_name << "initial C_" << i;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(InitialComposition,
                                        "initial composition",
                                        "Implementation of a plugin in which the tracer "
                                        "property is given as the initial composition "
                                        "at this position. The tracer gets as many "
                                        "properties as compositional fields."
                                        "\n\n")
    }
  }
}

