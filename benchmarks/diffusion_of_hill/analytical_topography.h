/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_analytical_topography_h
#define _aspect_postprocess_analytical_topography_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>


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
        virtual
        std::pair<std::string,std::string> execute (TableHandler &statistics);

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
        virtual
        void
        parse_parameters (ParameterHandler &prm);
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
        double last_output_time;

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
