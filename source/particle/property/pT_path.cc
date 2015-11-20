/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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
      PTPath<dim>::initialize_one_particle_property(const Point<dim> &,
                                                    const Vector<double> &solution,
                                                    const std::vector<Tensor<1,dim> > &,
                                                    std::vector<double> &data) const
      {
        data.push_back(solution[this->introspection().component_indices.pressure]);
        data.push_back(solution[this->introspection().component_indices.temperature]);
      }

      template <int dim>
      void
      PTPath<dim>::update_one_particle_property(const unsigned int data_position,
                                                const Point<dim> &,
                                                const Vector<double> &solution,
                                                const std::vector<Tensor<1,dim> > &,
                                                std::vector<double> &data) const
      {
        data[data_position] = solution[this->introspection().component_indices.pressure];
        data[data_position+1] = solution[this->introspection().component_indices.temperature];
      }

      template <int dim>
      UpdateTimeFlags
      PTPath<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      PTPath<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("p",1));
        property_information.push_back(std::make_pair("T",1));
        return property_information;
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
                                        "property is defined as the current pressure and "
                                        "temperature at this position. This can be used "
                                        "to generate pressure-temperature paths of "
                                        "material points over time."
                                        "\n\n")
    }
  }
}

