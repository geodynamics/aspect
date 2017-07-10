/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_adiabatic_conditions_function_h
#define _aspect_adiabatic_conditions_function_h


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
     * a given function with three components: temperature, pressure, density.
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
         * Initialization function.
         */
        virtual void initialize ();

        /**
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated, this is always true for the Function class.
         */
        virtual bool is_initialized() const;

        /**
         * Return the adiabatic temperature at a given point of the domain.
         */
        virtual double temperature (const Point<dim> &p) const;


        /**
         * Return the adiabatic pressure at a given point of the domain.
         */
        virtual double pressure (const Point<dim> &p) const;

        /**
         * Return the reference_density at a given point of the domain.
         */
        virtual
        double density (const Point<dim> &p) const;

        /**
         * Return the derivative of the density with respect to depth
         * at the given point @p p.
         */
        virtual
        double density_derivative (const Point<dim> &p) const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * ParsedFunction: depth->(temperature, pressure, density)
         */
        Functions::ParsedFunction<1> function;
    };
  }
}


#endif
