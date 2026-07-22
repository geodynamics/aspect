/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_plateness_statistics_h
#define _aspect_postprocess_plateness_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that computes surface plateness diagnostics on the top
     * boundary using the second invariant of the deviatoric strain-rate
     * tensor. F80 and F90 are the smallest fractions of the top-boundary area
     * that contain 80% and 90% of the integrated strain-rate invariant,
     * respectively. Consequently, smaller values indicate more strongly
     * localized surface deformation.
     *
     * The corresponding plateness values are computed as
     * p = 1 - F/reference_fraction. A value of one represents the limiting
     * case of deformation localized into an infinitesimally small area, zero
     * corresponds to the chosen reference fraction, and negative values
     * indicate deformation that is more distributed than the reference case.
     * Plateness is a relative diagnostic whose interpretation depends on the
     * reference fraction, model setup, and numerical resolution.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class PlatenessStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the surface plateness statistics and add them to the
         * statistics object.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

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

      private:
        /**
         * The reference surface-area fraction used to compute plateness as
         * p = 1 - F/reference_fraction. Changing this value changes the zero
         * point of p80 and p90, but not the computed F80 and F90 values.
         */
        double reference_fraction;
    };
  }
}


#endif
