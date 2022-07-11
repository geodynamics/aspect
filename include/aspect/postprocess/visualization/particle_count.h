/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_particle_count_h
#define _aspect_postprocess_visualization_particle_count_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from CellDataVectorCreator that takes an output
       * vector and computes a variable that represents the number of particles
       * in each cell.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class ParticleCount
        : public CellDataVectorCreator<dim>,
          public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           */
          ParticleCount();

          /**
           * @copydoc CellDataVectorCreator<dim>::execute()
           */
          std::pair<std::string, Vector<float> *>
          execute () const override;

          /**
           * Let the postprocessor manager know about the other postprocessors
           * this one depends on. Specifically, the particles postprocessor.
           */
          std::list<std::string>
          required_other_postprocessors() const override;
      };
    }
  }
}

#endif
