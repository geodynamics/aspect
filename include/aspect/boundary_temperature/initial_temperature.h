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


#ifndef __aspect__boundary_temperature_initial_temperature_h
#define __aspect__boundary_temperature_initial_temperature_h

#include <aspect/boundary_temperature/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace BoundaryTemperature
  {
    /**
     * A class that implements a temperature boundary condition for an
     * arbitrary geometry in which the temperature at the boundaries are the
     * same as in the initial conditions.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class InitialTemperature : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns the boundary temperatures that are defined
         * by the initial conditions.
         *
         * @copydoc aspect::BoundaryTemperature::Interface::boundary_temperature()
         */
        virtual
        double boundary_temperature (const types::boundary_id boundary_indicator,
                                     const Point<dim> &location) const;

        /**
         * Return the minimal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Return the maximal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary temperatures.
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
         * Temperatures at the inner and outer boundaries.
         */
        double min_temperature;
        double max_temperature;
    };
  }
}


#endif
