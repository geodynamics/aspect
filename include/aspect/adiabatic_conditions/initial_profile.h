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


#ifndef __aspect__adiabatic_conditions_initial_profile_h
#define __aspect__adiabatic_conditions_initial_profile_h


#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator.h>
#include <deal.II/base/point.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    using namespace dealii;

    /**
     * A class that represents adiabatic conditions, i.e. that starts at the top
     * of the domain and integrates pressure and temperature as we go down into
     * depth.
     *
     * @note The implementation has numerous deficiencies indicated in the .cc
     * file and may not quite compute what we want. Specifically, it doesn't
     * currently take into account all the physical parameters it needs, and it
     * also doesn't get gravity right with the exception of the simplest cases.
     */
    template <int dim>
    class InitialProfile : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
      /**
       * Constructor. Compute the adiabatic conditions along a vertical
       * transect of the geometry based on the given material model and other
       * quantities.
       */
      InitialProfile ();

      void update ();

      void initialize ();

      bool is_initialized() const;

      /**
       * Return the adiabatic temperature at a given point of the domain.
       */
      double temperature (const Point<dim> &p) const;

      /**
       * Return the adiabatic temperature profile as a vector of values
       * corresponding to increasing depth.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void get_adiabatic_temperature_profile(std::vector<double> &values) const;

      /**
       * Return the adiabatic pressure at a given point of the domain.
       */
      double pressure (const Point<dim> &p) const;

      private:

        /**
         * Wether the adiabatic conditions are already calculated.
         * This is important for plugins that are used by the adiabatic conditions
         * but also depend on the adiabatic conditions. This way they can behave
         * differently in initialization and model run.
         */
        bool initialized;

        /**
         * Number of points at which we compute the adiabatic values.
         */
        const unsigned int n_points;

        /**
         * Vectors of values of temperatures and pressures on a transect into
         * depth at which we have computed them. The public member functions of
         * this class interpolate linearly between these points.
         */
        std::vector<double> temperatures, pressures;

        /**
         * Interval spacing between each two data points in the tables above
         * with regard to the depth coordinate.
         */
        double delta_z;
    };
  }
}


#endif
