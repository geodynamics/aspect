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


#ifndef __aspect__initial_conditions_box_h
#define __aspect__initial_conditions_box_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes a perturbed initial temperature field for a box
     * geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class PerturbedBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;
    };

    /**
     * A class that describes an opposing poles initial temperature field for
     * a box geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class PolarBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;
    };

    /**
     * A field that describes a fractal initial temperature field
     *
     * @ingroup InitialCOnditionsModels
     */
    template <int dim>
    class MandelBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;
    };

    /**
     * A class that describes a shaped inclusion initial temperature field for
     * a box geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class InclusionShapeBox : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature(const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        std::string inclusion_shape;
        std::string inclusion_gradient;
        double radius;
        double ambient_temperature;
        double inclusion_temperature;
        double center_x;
        double center_y;
        double center_z;
    };
  }
}

#endif
