/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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
        void initialize () override;

        /**
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated, this is always true for the Function class.
         */
        bool is_initialized() const override;

        /**
         * Return the adiabatic temperature at a given point of the domain.
         */
        double temperature (const Point<dim> &p) const override;


        /**
         * Return the adiabatic pressure at a given point of the domain.
         */
        double pressure (const Point<dim> &p) const override;

        /**
         * Return the reference_density at a given point of the domain.
         */
        double density (const Point<dim> &p) const override;

        /**
         * Return the derivative of the density with respect to depth
         * at the given point @p p.
         */
        double density_derivative (const Point<dim> &p) const override;

        /**
         * Declare the parameters for the input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:

        /**
         * ParsedFunction: depth->(temperature, pressure, density)
         */
        Functions::ParsedFunction<1> function;
    };
  }
}


#endif
