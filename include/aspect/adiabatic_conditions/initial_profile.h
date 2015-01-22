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
     * A simple class that calculates adiabatic conditions. This
     * implementation calculates a simple profile at model start time and does
     * not update it over time. It utilizes the initial condition for
     * compositional fields at a reference point (a generic point in the
     * current depth) to calculate the material parameters. The gravity is
     * assumed to be directly downward, i.e. radial in spherical models and
     * vertical in box-shaped models.
     */
    template <int dim>
    class InitialProfile : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        InitialProfile ();

        /**
         * Initialization function. Because this function is called after
         * initializing the SimulatorAccess, all of the necessary information
         * is available to calculate the adiabatic profile. It computes the
         * adiabatic conditions along a vertical transect of the geometry
         * based on the given material model and other quantities.
         */
        virtual void initialize ();

        /**
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated. This is for example the case for the simple
         * com- pressible material model, which uses the adiabatic temperature
         * as reference temperature to calculate the density. For the
         * calculation of the adiabatic conditions this functionality is
         * simply switched off, because we are always on the reference
         * profile. This way the plugin behaves differently at initialization
         * time of the adiabatic conditions and during the main model run.
         */
        virtual bool is_initialized() const;

        /**
         * Empty update function. This class does not update the adiabatic
         * profile over time.
         */
        virtual void update ();

        /**
         * Return the adiabatic temperature at a given point of the domain.
         */
        virtual double temperature (const Point<dim> &p) const;

        /**
         * Return the adiabatic temperature profile as a vector of values
         * corresponding to increasing depth.
         *
         * @param values The output vector of depth averaged values. The
         * function takes the pre-existing size of this vector as the number
         * of depth slices.
         */
        virtual void get_adiabatic_temperature_profile(std::vector<double> &values) const;

        /**
         * Return the adiabatic pressure at a given point of the domain.
         */
        virtual double pressure (const Point<dim> &p) const;

      private:

        /**
         * Wether the adiabatic conditions are already calculated. This is
         * important for plugins that are used by the adiabatic conditions but
         * also depend on the adiabatic conditions. This way they can behave
         * differently in initialization and model run.
         */
        bool initialized;

        /**
         * Number of points at which we compute the adiabatic values.
         */
        const unsigned int n_points;

        /**
         * Vectors of values of temperatures and pressures on a transect into
         * depth at which we have computed them. The public member functions
         * of this class interpolate linearly between these points.
         */
        std::vector<double> temperatures;
        std::vector<double> pressures;

        /**
         * Interval spacing between each two data points in the tables above
         * with regard to the depth coordinate.
         */
        double delta_z;
    };
  }
}


#endif
