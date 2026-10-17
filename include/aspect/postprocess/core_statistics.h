// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_core_statistics_h
#define _aspect_postprocess_core_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the dynamic core.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class CoreStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Evaluate statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

      private:
        /**
         * Controls whether output the total excess entropy or the individual entropy terms
         * (i.e. entropy for specific heat, radioactive heating, gravitational contribution,
         * and adiabatic contribution).
         *
         * This variable is read from the parameter file through a parameter called 'Excess entropy only'.
         */
        bool   excess_entropy_only;
    };
  }
}


#endif
