/*
  Copyright (C) 2014 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_boundary_temperature_constant_h
#define _aspect_boundary_temperature_constant_h

#include <aspect/boundary_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <map>


namespace aspect
{
  namespace BoundaryTemperature
  {
    /**
     * A class that implements a temperature boundary condition for an
     * arbitrary geometry, where it returns a constant value for a given
     * boundary indicator. Takes an ordered list of temperatures, where the
     * position in the list corresponds to the indicator for that temperature.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class Constant : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns user-defined constant temperatures at the
         * boundaries of arbitrary geometries.
         *
         * @copydoc aspect::BoundaryTemperature::Interface::boundary_temperature()
         */
        double boundary_temperature (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position) const override;

        /**
         * Return the minimal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const override;

        /**
         * Return the maximal the temperature on that part of the boundary on
         * which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const override;

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
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Temperatures at the boundaries.
         */
        std::map<types::boundary_id, double> boundary_temperatures;
    };
  }
}


#endif
