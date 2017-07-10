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
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */


#include <aspect/particle/property/composition.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      Composition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                         std::vector<double> &data) const
      {
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          data.push_back(this->get_initial_composition_manager().initial_composition(position,i));
      }

      template <int dim>
      void
      Composition<dim>::update_one_particle_property(const unsigned int data_position,
                                                     const Point<dim> &,
                                                     const Vector<double> &solution,
                                                     const std::vector<Tensor<1,dim> > &,
                                                     const ArrayView<double> &data) const
      {
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            const unsigned int solution_component = this->introspection().component_indices.compositional_fields[i];
            data[data_position+i] = solution[solution_component];
          }
      }

      template <int dim>
      UpdateTimeFlags
      Composition<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      Composition<dim>::get_needed_update_flags () const
      {
        return update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      Composition<dim>::get_property_information() const
      {

        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <composition>, "
                               "but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));

        std::vector<std::pair<std::string,unsigned int> > property_information;



        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            const std::string field_name = this->introspection().name_for_compositional_index(i);
            property_information.push_back(std::make_pair(field_name,1));
          }
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(Composition,
                                        "composition",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined by the compositional fields in "
                                        "the model. This can be used to track solid composition"
                                        "evolution over time.")
    }
  }
}

