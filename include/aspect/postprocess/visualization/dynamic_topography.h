/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_dynamic_topography_h
#define _aspect_postprocess_visualization_dynamic_topography_h

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
       * vector and computes a variable that represents the dynamic
       * topography. This quantity, strictly speaking, only makes sense at the
       * surface of the domain. Thus, the value is set to zero in all the
       * cells inside of the domain.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class DynamicTopography
        : public CellDataVectorCreator<dim>,
          public SimulatorAccess<dim>
      {
        public:
          /**
           * @copydoc CellDataVectorCreator<dim>::execute()
           */
          virtual
          std::pair<std::string, Vector<float> *>
          execute () const;

          /**
           * Register the other postprocessor that we need: DynamicTopography
           */
          virtual
          std::list<std::string>
          required_other_postprocessors() const;
      };
    }
  }
}

#endif
