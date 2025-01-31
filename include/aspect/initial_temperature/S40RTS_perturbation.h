/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_S40RTS_perturbation_h
#define _aspect_initial_temperature_S40RTS_perturbation_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialTemperature
  {
    namespace internal
    {
      namespace S40RTS
      {
        class SphericalHarmonicsLookup
        {
          public:
            SphericalHarmonicsLookup(const std::string &filename,
                                     const MPI_Comm comm);

            /// Declare a function that returns the cosine coefficients
            const std::vector<double> &
            cos_coeffs() const;

            /// Declare a function that returns the sine coefficients
            const std::vector<double> &
            sin_coeffs() const;

            unsigned int maxdegree() const;

          private:
            unsigned int order;
            std::vector<double> a_lm;
            std::vector<double> b_lm;
        };

        class SplineDepthsLookup
        {
          public:
            SplineDepthsLookup(const std::string &filename,
                               const MPI_Comm comm);

            const std::vector<double> &
            spline_depths() const;

          private:
            std::vector<double> depths;
        };
      }
    }

    template <int dim>
    class PatchOnS40RTS;


    /**
     * A class that describes a perturbed initial temperature field for a
     * spherical shell geometry model. The perturbation is based on the S20RTS
     * / S40RTS global shear wave velocity model by Ritsema et al.
     * http://www.earth.lsa.umich.edu/~jritsema/research.html
     *
     * @ingroup InitialTemperatures
     */

    template <int dim>
    class S40RTSPerturbation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        S40RTSPerturbation ();

        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        void
        initialize () override;

        /**
         * Return the Vs as a function of position.
         */
        virtual
        double get_Vs (const Point<dim> &position) const;

        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

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
         * An enum to describe which method should be chosen to scale vs to density.
         */
        enum VsToDensityMethod
        {
          file,
          constant
        };

        /**
         * Currently chosen source for vs to density scaling.
         */
        VsToDensityMethod vs_to_density_method;

        /**
         * File directory and names
         */
        std::string data_directory;
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
         * velocity into a temperatures perturbation. The first parameter is
         * constant so far but could be made depth dependent as constraint by
         * e.g. Forte, A.M. & Woodward, R.L., 1997. Seismic-geodynamic
         * constraints on three- dimensional structure, vertical flow, and
         * heat transfer in the mantle, J. Geophys. Res. 102 (B8),
         * 17,981-17,994.
         * The last parameter is a depth down to which heterogeneities are
         * zeroed out.
         */
        double vs_to_density_constant;
        double thermal_alpha;
        double no_perturbation_depth;

        /**
         * This parameter allows to remove the degree 0 component of the shear
         * wave velocity perturbation, which guarantees that average
         * temperature at a certain depth is the background temperature.
         */
        bool zero_out_degree_0;

        /**
         * This parameter allows to use a lower maximum degree when reading
         * the spherical harmonic data file.
         */
        bool lower_max_degree;

        /**
         * The maximum degree the users specify, which is only valid when
         * "lower_max_degree" is set to true.
         */
        unsigned int specified_max_degree;

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
        std::unique_ptr<internal::S40RTS::SphericalHarmonicsLookup> spherical_harmonics_lookup;

        /**
         * Pointer to an object that reads and processes the depths for the
         * spline knot points.
         */
        std::unique_ptr<internal::S40RTS::SplineDepthsLookup> spline_depths_lookup;

        /**
         * Object containing the data profile.
         */
        aspect::Utilities::AsciiDataProfile<dim> profile;

        /**
         * The column index of the vs to density scaling in the data file
         */
        unsigned int vs_to_density_index;

        /**
         * Whether to use the thermal expansion coefficient from the material model
         */
        bool use_material_model_thermal_alpha;

        template <int dim2> friend class PatchOnS40RTS;
    };

  }
}

#endif
