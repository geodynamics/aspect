// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}


#endif
