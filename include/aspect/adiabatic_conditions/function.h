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


#ifndef __aspect__adiabatic_conditions_function_h
#define __aspect__adiabatic_conditions_function_h


#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parsed_function.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    using namespace dealii;

    /**
     * A simple class that sets the adiabatic conditions based on given
     * functions.
     */
    template <int dim>
    class Function : public Interface<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        Function ();

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
         * Return the adiabatic pressure at a given point of the domain.
         */
        virtual double pressure (const Point<dim> &p) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * ParsedFunction: depth->(temperature, pressure)
         */
        Functions::ParsedFunction<1> function;
    };
  }
}


#endif
