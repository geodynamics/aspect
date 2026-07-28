/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_box_h
#define _aspect_initial_temperature_box_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace InitialTemperature
  {
    /**
     * A class that describes a perturbed initial temperature field for a box
     * geometry.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class PerturbedBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;
    };

    /**
     * A class that describes an opposing poles initial temperature field for
     * a box geometry.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class PolarBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;
    };

    /**
     * A field that describes a fractal initial temperature field
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class MandelBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;
    };

    /**
     * A class that describes a shaped inclusion initial temperature field for
     * a box geometry.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class InclusionShapeBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature(const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * This variable is read from the parameter file through a parameter called 'Inclusion shape'.
         */
        std::string inclusion_shape;
        /**
         * This variable is read from the parameter file through a parameter called 'Inclusion gradient'.
         */
        std::string inclusion_gradient;
        /**
         * This variable is read from the parameter file through a parameter called 'Shape radius'.
         */
        double radius;
        /**
         * This variable is read from the parameter file through a parameter called 'Ambient temperature'.
         */
        double ambient_temperature;
        /**
         * This variable is read from the parameter file through a parameter called 'Inclusion temperature'.
         */
        double inclusion_temperature;
        /**
         * This variable is read from the parameter file through a parameter called 'Center X'.
         */
        double center_x;
        /**
         * This variable is read from the parameter file through a parameter called 'Center Y'.
         */
        double center_y;
        /**
         * This variable is read from the parameter file through a parameter called 'Center Z'.
         */
        double center_z;
    };
  }
}

#endif
