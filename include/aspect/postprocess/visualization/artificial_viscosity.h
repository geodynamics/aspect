/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_visualization_artificial_viscosity_h
#define __aspect__postprocess_visualization_artificial_viscosity_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived that implements a function that provides the
       * artificial viscosity for the temperature equation on each cell for
       * graphical output.
       */
      template <int dim>
      class ArtificialViscosity : public CellDataVectorCreator<dim>,
        public SimulatorAccess<dim>
      {
        public:
          /**
           * The function classes have to implement that want to output
           * cellwise data.
           * @return A pair of values with the following meaning: - The first
           * element provides the name by which this data should be written to
           * the output file. - The second element is a pointer to a vector
           * with one element per active cell on the current processor.
           * Elements corresponding to active cells that are either artificial
           * or ghost cells (in deal.II language, see the deal.II glossary)
           * will be ignored but must nevertheless exist in the returned
           * vector. While implementations of this function must create this
           * vector, ownership is taken over by the caller of this function
           * and the caller will take care of destroying the vector pointed
           * to.
           */
          virtual
          std::pair<std::string, Vector<float> *>
          execute () const;
      };
    }
  }
}

#endif
