/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_finite_element_information_h
#define _aspect_postprocess_finite_element_information_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that prints the finite element space used by every
     * solution variable at time zero.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class FiniteElementInformation : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Print the finite element spaces used by the solution variables.
         */
        std::pair<std::string, std::string>
        execute(TableHandler &statistics) override;
    };
  }
}


#endif
