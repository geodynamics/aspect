/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__boundary_temperature_constant_h
#define __aspect__boundary_temperature_constant_h

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
         * Return the temperature that is to hold at a particular location on
         * the boundary of the domain. This function returns the constant
         * temperatures read from the parameter file.
         *
         * @param geometry_model The geometry model that describes the domain.
         * This may be used to determine whether the boundary temperature
         * model is implemented for this geometry.
         * @param boundary_indicator The boundary indicator of the part of the
         * boundary of the domain on which the point is located at which we
         * are requesting the temperature.
         * @param location The location of the point at which we ask for the
         * temperature.
         */
        virtual
        double temperature (const GeometryModel::Interface<dim> &geometry_model,
                            const unsigned int                   boundary_indicator,
                            const Point<dim>                    &location) const;

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
        std::map<types::boundary_id, double> boundary_temperatures;
    };
  }
}


#endif
