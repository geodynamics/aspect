/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#include <aspect/particle/property/melt_particle.h>
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
      Melt<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                  std::vector<double> &data) const
      {
        data.push_back(0.0);
      }

      template <int dim>
      void
      Melt<dim>::update_one_particle_property(const unsigned int data_position,
                                              const Point<dim> &,
                                              const Vector<double> &solution,
                                              const std::vector<Tensor<1,dim> > &,
                                              const ArrayView<double> &data) const
      {
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("Material model Melt simple with melt transport only "
                               "works if there is a compositional field called porosity."));
        const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
        const unsigned int solution_component = this->introspection().component_indices.compositional_fields[porosity_idx];

        if (solution[this->introspection().component_indices.compositional_fields[porosity_idx]] > this->get_melt_handler().melt_transport_threshold)
          data[data_position] = 1.0;
        else
          data[data_position] =  0.0;
      }

      template <int dim>
      UpdateTimeFlags
      Melt<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      Melt<dim>::get_needed_update_flags () const
      {
        return update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      Melt<dim>::get_property_information() const
      {
        std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("porosity",1));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(Melt,
                                        "melt",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as presence of melt above the "
                                        "melt transport theshold. Returns a porosity of 0 "
                                        "if melt is not present and a porosity of 1 if melt  "
                                        "is present.")
    }
  }
}

