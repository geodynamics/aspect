/*
  Copyright (C) 2016 - 2018 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/adiabatic_boundary.h>
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  namespace InitialTemperature
  {


    template <int dim>
    AdiabaticBoundary<dim>::AdiabaticBoundary ()
      :
      surface_boundary_id(1)
    {}

    template <int dim>
    void
    AdiabaticBoundary<dim>::initialize ()
    {
      // Find the boundary indicator that represents the surface
      surface_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("outer");

      std::set<types::boundary_id> surface_boundary_set;
      surface_boundary_set.insert(surface_boundary_id);

      // The input ascii table contains two components, the crust depth and the LAB depth
      Utilities::AsciiDataBoundary<dim>::initialize(surface_boundary_set,
                                                    2);
    }

    template <int dim>
    double
    AdiabaticBoundary<dim>::initial_temperature (const Point<dim> &position) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double isotherm_depth              =
        Utilities::AsciiDataBoundary<dim>::get_data_component(surface_boundary_id, position, 0);
      //std::cout<<isotherm_depth<<std::endl;
      if (depth > isotherm_depth)
        return isotherm_temperature + (depth - isotherm_depth) * temperature_gradient;
      else
        return surface_temperature + (depth/isotherm_depth) * (isotherm_temperature - surface_temperature);
    }

    template <int dim>
    void
    AdiabaticBoundary<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/",
                                                          "litho.kenya.txt");
        prm.enter_subsection("Adiabatic boundary");
        {
          prm.declare_entry ("Isotherm temperature", "1673.15",
                             Patterns::Double (0),
                             "The value of the isothermal boundary temperature. Units: Kelvin.");
          prm.declare_entry ("Surface temperature", "273.15",
                             Patterns::Double (0),
                             "The value of the surface temperature. Units: Kelvin.");
          prm.declare_entry ("Adiabatic temperature gradient", "0.0005",
                             Patterns::Double (0),
                             "The value of the adiabatic temperature gradient. Units: $K m^{-1}$.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    AdiabaticBoundary<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        Utilities::AsciiDataBase<dim>::parse_parameters(prm);

        prm.enter_subsection("Adiabatic boundary");
        {
          isotherm_temperature = prm.get_double("Isotherm temperature");
          surface_temperature  = prm.get_double("Surface temperature");
          temperature_gradient = prm.get_double("Adiabatic temperature gradient");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }

  }

}


namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(AdiabaticBoundary,
                                              "adiabatic boundary",
                                              "An initial temperature condition that allows for discretizing "
                                              "the model domain into two layers separated by a user-defined "
                                              "isothermal boundary using a table look-up approach. The user includes an "
                                              "input ascii data file that is formatted as 3 columns of `latitude', "
                                              "`longitude', and `depth', where `depth' is in kilometers and "
                                              "represents the depth of an initial temperature of 1673.15 K (by default). "
                                              "The temperature is defined from the surface (273.15 K) to the isotherm "
                                              "as a linear gradient. Below the isotherm the temperature increases "
                                              "approximately adiabatically (0.0005 K per meter). This initial temperature condition "
                                              "is designed specifically for the ellipsoidal chunk geometry model.")
  }
}
