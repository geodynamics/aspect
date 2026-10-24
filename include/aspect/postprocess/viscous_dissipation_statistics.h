// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2021 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_viscous_dissipation_statistics_h
#define _aspect_postprocess_viscous_dissipation_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes the viscous dissipation over the area
     * spanned by each compositional field and over the whole domain.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class ViscousDissipationStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution to compute the viscous dissipation
         * per compositional field and over the whole domain.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
    };
  }
}


#endif
