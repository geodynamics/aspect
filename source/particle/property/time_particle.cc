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


#include <aspect/particle/property/time_particle.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
	
      template <int dim>
      void
      TimeParticle<dim>::initialize_one_particle_property(const Point<dim> &/*position*/,
                                                  std::vector<double> &data) const
													  {
	  data.push_back(0.0);
													  }

      template <int dim>
      void
      TimeParticle<dim>::update_one_particle_property(const unsigned int data_position,
                                                  const Point<dim> &,
                                                  const Vector<double> &time,
                                                  const std::vector<Tensor<1,dim> > &,
                                                  const ArrayView<double> &data) const
												  {
      data[data_position] = this->get_time();
												  }
      
      template <int dim>
      UpdateTimeFlags
      TimeParticle<dim>::need_update() const
      {
        return update_output_step;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int> >
      TimeParticle<dim>::get_property_information() const
      {
        /*const */std::vector<std::pair<std::string,unsigned int> > property_information (1,std::make_pair("time",1));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(TimeParticle,
                                        "time",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined as the current time.")
    }
  }
}

