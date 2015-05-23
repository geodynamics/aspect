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


#ifndef __aspect__postprocess_geoid_h
#define __aspect__postprocess_geoid_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    using namespace dealii;

    namespace internal
    {
      /**
       * A struct that contains the representation of the coefficients of a
       * spherical harmonic expansion.
       */
      struct HarmonicCoefficients
      {
          HarmonicCoefficients(const unsigned int max_degree);

          std::vector<double> sine_coefficients;
          std::vector<double> cosine_coefficients;
      };

      /**
       * A class to expand arbitrary fields of doubles into spherical harmonic
       * coefficients
       */
      template <int dim>
      class SphericalHarmonicsExpansion
      {
        public:
          SphericalHarmonicsExpansion(const unsigned int max_degree);

          void add_data_point (const Point<dim> &position,
                         const double value);

          HarmonicCoefficients
          get_coefficients () const;

          void
          mpi_sum_coefficients (MPI_Comm mpi_communicator);

        private:
          const unsigned int max_degree;

          HarmonicCoefficients coefficients;
      };
    }

    /**
     * A postprocessor that computes the geoid height at the surface.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Geoid : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for the geoid height at the surface.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);


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
         * A parameter that we read from the input file that denotes whether
         * we should include the contribution by dynamic topography.
         */
        bool include_topography_contribution;

        double core_mass;
        double density_below;
        double density_above;
        unsigned int max_degree;

        /**
         * The geoid contribution is added on a per-layer basis. These are the
         * coefficients for each layer, which will be finally added to a
         * combined contribution at the surface.
         */
        std::vector <std_cxx11::shared_ptr<internal::SphericalHarmonicsExpansion<dim> > > internal_density_expansion;
        std_cxx11::shared_ptr<internal::SphericalHarmonicsExpansion<dim> > surface_topography_expansion;
        std_cxx11::shared_ptr<internal::SphericalHarmonicsExpansion<dim> > bottom_topography_expansion;
    };
  }
}


#endif
