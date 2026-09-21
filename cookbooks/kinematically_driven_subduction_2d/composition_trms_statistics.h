// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2024 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_composition_trms_statistics_h
#define _aspect_postprocess_composition_trms_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes the root mean square temperature
     * over the area spanned by each of the compositional fields.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class CompositionTrmsStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution to compute the root-mean-square temperature per compositional field.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static void
        declare_parameters(ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters(ParameterHandler &prm) override;
        /**
         * @}
         */

      private:
        std::vector<std::string> selected_fields;
    };
  }
}


#endif
