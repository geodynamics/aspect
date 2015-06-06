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

#include <aspect/particle/output/null.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      template <int dim>
      NullOutput<dim>::NullOutput()
        :
        Interface<dim> ()
      {}

      template <int dim>
      std::string
      NullOutput<dim>::output_particle_data(const std::multimap<LevelInd, Particle<dim> > &/*particles*/,
                                            const std::vector<std::string>  &/*data_names*/,
                                            const std::vector<unsigned int> &/*data_components*/,
                                            const double &/*current_time*/)
      {
        return "";
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      ASPECT_REGISTER_PARTICLE_OUTPUT(NullOutput,
                                      "none",
                                      "Blank class to avoid output of data, for example, if "
                                      "particles are used to represent other physical phenomenon "
                                      "that influences the simulation and we don't care about "
                                      "their positions.")
    }
  }
}


