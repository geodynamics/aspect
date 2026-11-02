// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2020-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_analytical_topography_h
#define _aspect_postprocess_analytical_topography_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <limits>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs the surface topography to file.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class AnalyticalTopography : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Output predicted and analytical topography [m] to file.
         */
        std::pair<std::string,std::string> execute (TableHandler &statistics) override;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
         * @}
         */

      private:
        /**
         * Whether or not to produce text files with topography values
         */
        bool write_to_file;

        /**
         * Interval between the generation of text output. This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last text output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        double last_output_time = std::numeric_limits<double>::lowest();

        /**
         * The amplitude of the sinusoidal initial topography.
         */
        double amplitude;

        /**
         * The diffusivity used in computing the diffusion of the surface
         * topography.
         */
        double kappa;

        /**
         * The width (X extent) of the 2D rectangular domain.
         */
        double domain_width;

        /**
         * A switch to vary between different analytical solutions.
         */
        unsigned int analytical_solution_example;

    };
  }
}


#endif
