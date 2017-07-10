/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#ifndef _aspect__adiabatic_conditions_ascii_data_h
#define _aspect__adiabatic_conditions_ascii_data_h


#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/point.h>


namespace aspect
{
  namespace AdiabaticConditions
  {
    using namespace dealii;

    /**
     * A simple class that reads adiabatic conditions from a file.
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataProfile<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        AsciiData ();

        /**
         * Initialization function. Because this function is called after
         * initializing the SimulatorAccess, all of the necessary information
         * is available to calculate the adiabatic profile. It computes the
         * adiabatic conditions along a vertical transect of the geometry
         * based on the given material model and other quantities.
         */
        virtual void initialize ();

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataProfile<dim>::initialize;

        /**
         * Some plugins need to know whether the adiabatic conditions are
         * already calculated. This is for example the case for the simple
         * compressible material model, which uses the adiabatic temperature
         * as reference temperature to calculate the density. For the
         * calculation of the adiabatic conditions this functionality is
         * simply switched off, because we are always on the reference
         * profile. This way the plugin behaves differently at initialization
         * time of the adiabatic conditions and during the main model run.
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
         * Return the reference density at a given point of the domain.
         */
        virtual double density (const Point<dim> &p) const;

        /**
         * Return the derivative of the density with respect to depth
         * at the given point @p p.
         */
        virtual
        double density_derivative (const Point<dim> &p) const;


        /**
         * Declare the parameters for the input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * Whether the adiabatic conditions are already calculated. This is
         * important for plugins that are used by the adiabatic conditions but
         * also depend on the adiabatic conditions. This way they can behave
         * differently in initialization and model run.
         */
        bool initialized;

        /**
         * The column indices of the temperature, pressure, and density column
         * in the data file.
         */
        unsigned int temperature_index;
        unsigned int pressure_index;
        unsigned int density_index;
    };
  }
}


#endif
