/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#include <deal.II/base/signaling_nan.h>
#include <aspect/boundary_temperature/spherical_constant.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>

#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryTemperature
  {
// ------------------------------ SphericalConstant -------------------

    template <int dim>
    double
    SphericalConstant<dim>::
    boundary_temperature (const types::boundary_id boundary_indicator,
                          const Point<dim> &) const
    {
      Assert (this->get_geometry_model().translate_id_to_symbol_name(boundary_indicator) == "top" ||
              this->get_geometry_model().translate_id_to_symbol_name(boundary_indicator) == "bottom",
              ExcMessage ("Unknown boundary indicator for geometry model. "
                          "The given boundary should be ``top'' or ``bottom''."));

      return boundary_temperatures[boundary_indicator];
    }


    template <int dim>
    double
    SphericalConstant<dim>::
    minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      double min = std::numeric_limits<double>::max();

      for (unsigned int id=0; id<boundary_temperatures.size(); ++id)
        if (fixed_boundary_ids.empty() || fixed_boundary_ids.find(id) != fixed_boundary_ids.end())
          if (std::isnan(boundary_temperatures[id]) == false)
            min = std::min(min,boundary_temperatures[id]);

      return min;
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const
    {
      double max = -std::numeric_limits<double>::max();

      for (unsigned int id=0; id<boundary_temperatures.size(); ++id)
        if (std::isnan(boundary_temperatures[id]) == false)
          if (fixed_boundary_ids.empty() || fixed_boundary_ids.find(id) != fixed_boundary_ids.end())
            max = std::max(max,boundary_temperatures[id]);

      return max;
    }



    template <int dim>
    void
    SphericalConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Spherical constant");
        {
          prm.declare_entry ("Outer temperature", "0.",
                             Patterns::Double (),
                             "Temperature at the outer boundary (lithosphere water/air). Units: \\si{\\kelvin}.");
          prm.declare_entry ("Inner temperature", "6000.",
                             Patterns::Double (),
                             "Temperature at the inner boundary (core mantle boundary). Units: \\si{\\kelvin}.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    SphericalConstant<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Spherical constant");
        {
          const auto boundary_names_map = this->get_geometry_model().get_symbolic_boundary_names_map();
          boundary_temperatures.resize(boundary_names_map.size(), std::numeric_limits<double>::quiet_NaN());

          if (boundary_names_map.find("bottom") != boundary_names_map.end())
            boundary_temperatures[boundary_names_map.at("bottom")] =  prm.get_double ("Inner temperature");

          if (boundary_names_map.find("top") != boundary_names_map.end())
            boundary_temperatures[boundary_names_map.at("top")] =  prm.get_double ("Outer temperature");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryTemperature
  {
    ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(SphericalConstant,
                                               "spherical constant",
                                               "A model in which the temperature is chosen constant on "
                                               "the inner and outer boundaries of a spherical shell, ellipsoidal "
                                               "chunk or chunk. "
                                               "Parameters are read from subsection 'Spherical constant'.")
  }
}
