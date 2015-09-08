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


#ifndef __aspect__initial_conditions_SAVANI_perturbation_h
#define __aspect__initial_conditions_SAVANI_perturbation_h

#include <aspect/simulator_access.h>
#include <deal.II/base/std_cxx11/array.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    namespace internal
    {
      namespace SAVANI
      {
        class SphericalHarmonicsLookup;
        class SplineDepthsLookup;
      }
    }

    /**
     * A class that describes a perturbed initial temperature field for a
     * spherical shell geometry model. The perturbation is based on the SAVANI
     * global shear wave velocity model by Auer et al.
     * http://n.ethz.ch/~auerl/research.html
     *
     * @ingroup InitialConditionsModels
     */

    template <int dim>
    class SAVANIPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        void
        initialize ();

        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;

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


      private:

        /**
         * File directory and names
         */
        std::string datadirectory;
        std::string spline_depth_file_name;

        /**
         * This parameter allows setting the input file for the SAVANI global
         * shear-wave perturbation.
         */
        std::string harmonics_coeffs_file_name;

        /**
         * The parameters below describe the perturbation of shear wave
         * velocity into a temperatures perturbation. The first parameter is
         * constant so far but could be made depth dependent as constraint by
         * e.g. Forte, A.M. & Woodward, R.L., 1997. Seismic-geodynamic
         * constraints on three- dimensional structure, vertical flow, and
         * heat transfer in the mantle, J. Geophys. Res. 102 (B8),
         * 17,981-17,994.
         * The last parameter is a depth down to which heterogeneities are
         * zeroed out.
         */
        double vs_to_density;
        double thermal_alpha;
        double no_perturbation_depth;

        /**
         * This parameter allows to remove the degree 0 component of the shear
         * wave velocity perturbation, which guarantees that average
         * temperature at a certain depth is the background temperature.
         */
        bool zero_out_degree_0;

        /**
         * This parameter gives the reference temperature, which will be
         * perturbed. In the compressional case the background temperature
         * will be the adiabat.
         */
        double reference_temperature;

        /**
         * Pointer to an object that reads and processes the spherical
         * harmonics coefficients
         */
        std_cxx11::shared_ptr<internal::SAVANI::SphericalHarmonicsLookup> spherical_harmonics_lookup;

        /**
         * Pointer to an object that reads and processes the depths for the
         * spline knot points.
         */
        std_cxx11::shared_ptr<internal::SAVANI::SplineDepthsLookup> spline_depths_lookup;

    };

  }
}

#endif
