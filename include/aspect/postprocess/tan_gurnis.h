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
/*  $Id$  */


#ifndef __aspect__postprocess_tan_gurnis_h
#define __aspect__postprocess_tan_gurnis_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that evaluates the accuracy of the solution of the
     * aspect::MaterialModel::TanGurnis material model.
     *
     * The implementation writes out the solution to be read in by a matlab
     * script. The benchmark is from the article
     * @code
     *  @article{tan2007compressible,
     *    title={Compressible thermochemical convection and application to lower mantle structures},
     *    author={Tan, E. and Gurnis, M.},
     *    journal={JOURNAL OF GEOPHYSICAL RESEARCH-ALL SERIES-},
     *    volume={112},
     *    number={B6},
     *    pages={6304},
     *    year={2007},
     *    publisher={AGU AMERICAN GEOPHYSICAL UNION}
     *    }
     * @endcode
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class TanGurnis : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
