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


#ifndef __aspect__initial_conditions_S40RTS_perturbation_h
#define __aspect__initial_conditions_S40RTS_perturbation_h

#include <aspect/simulator_access.h>
#include <deal.II/base/std_cxx1x/array.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    namespace internal
    {

      class SphericalHarmonicsLookup;
      class SplineDepthsLookup;
    }

    /**
     * A class that describes a perturbed initial temperature field for a
     * spherical shell geometry model. The perturbation is based on the S20RTS
     * / S40RTS global shear wave velocity model by Ritsema et al.
     * http://www.earth.lsa.umich.edu/~jritsema/research.html
     *
     * @ingroup InitialConditionsModels
     */

    template <int dim>
    class S40RTSPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * This parameter allows setting the input file for the shear-wave
         * perturbation. Options so far are S20RTS.sph and S40RTS.sph. For
         * S40RTS there are different versions available that differ by the
         * degree of damping in the seismic inversion. These models could be
         * downloaded and used as well.
         */
        std::string harmonics_coeffs_file_name;

        /**
         * The parameters below describe the perturbation of shear wave
         * velocity into a temperatures perturbation The first parameter is
         * constant so far but could be made depth dependent as constraint by
         * e.g. Forte, A.M. & Woodward, R.L., 1997. Seismic-geodynamic
         * constraints on three- dimensional structure, vertical flow, and
         * heat transfer in the mantle, J. Geophys. Res. 102 (B8),
         * 17,981-17,994.
         */
        double vs_to_density;
        double thermal_alpha;

        /**
         * This parameter allows to set the degree 0 component of the shear
         * wave velocity perturbation to zero, which guarantees that average
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
        std_cxx1x::shared_ptr<internal::SphericalHarmonicsLookup> spherical_harmonics_lookup;

        /**
         * Pointer to an object that reads and processes the depths for the
         * spline knot points.
         */
        std_cxx1x::shared_ptr<internal::SplineDepthsLookup> spline_depths_lookup;

    };

  }
}

#endif
