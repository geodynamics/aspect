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


#ifndef __aspect__postprocess_duretz_et_al_h
#define __aspect__postprocess_duretz_et_al_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that evaluates the accuracy of the solution of the
     * aspect::MaterialModel::DuretzEtAl::* material models.
     *
     * The implementation of error evaluators that correspond to the
     * benchmarks defined in the following paper:
     * @code
     *  @Article{DMGT11,
     *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
     *    title =        {Discretization errors and free surface stabilization in the
     *                  finite difference and marker-in-cell method for applied
     *                  geodynamics: {A} numerical study},
     *    journal =      {Geochemistry Geophysics Geosystems},
     *    year =         2011,
     *    volume =       12,
     *    pages =        {Q07004/1--26}}
     * @endcode
     *
     * @note While this paper summarizes the benchmarks used here, some of the
     * benchmarks actually originate in earlier papers. For the original
     * references, see the bibliography of the paper above.
     * @ingroup Postprocessing
     */
    template <int dim>
    class DuretzEtAl : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
