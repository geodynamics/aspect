/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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


#include <aspect/particle/property/peridotite.h>
#include <aspect/simulator.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      Peridotite<dim>::initialize_one_particle_property(const Point<dim> &/*position*/,
                                                        std::vector<double> &data) const
      {
        data.push_back(0.0);
      }

      template <int dim>
      void
      Peridotite<dim>::update_one_particle_property(const unsigned int data_position,
                                                    const Point<dim> &,
                                                    const Vector<double> &solution,
                                                    const std::vector<Tensor<1,dim> > &,
                                                    const ArrayView<double> &data) const
      {
        AssertThrow(this->introspection().compositional_name_exists("peridotite"),
                    ExcMessage("Particle property peridotite only works if "
                               "there is a compositional field called peridotite."));
        const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
        const unsigned int solution_component = this->introspection().component_indices.compositional_fields[peridotite_idx];

        data[data_position] = solution[this->introspection().component_indices.compositional_fields[peridotite_idx]];
      }

      template <int dim>
      UpdateTimeFlags
      Peridotite<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      Peridotite<dim>::get_needed_update_flags () const
      {
        return update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      Peridotite<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("peridotite",1));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(Peridotite,
                                        "peridotite",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the peridotite depletion "
                                        "at this position. This can be used to track solid "
                                        "composition evolution over time.")
    }
  }
}

