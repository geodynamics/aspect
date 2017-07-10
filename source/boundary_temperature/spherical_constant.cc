/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
      const GeometryModel::Interface<dim> *geometry_model = &this->get_geometry_model();
      const std::string boundary_name = geometry_model->translate_id_to_symbol_name(boundary_indicator);

      if (boundary_name == "bottom")
        return inner_temperature;
      else if (boundary_name =="top")
        return outer_temperature;
      else
        {
          Assert (false, ExcMessage ("Unknown boundary indicator for geometry model. "
                                     "The given boundary should be ``top'' or ``bottom''."));
          return numbers::signaling_nan<double>();
        }
    }


    template <int dim>
    double
    SphericalConstant<dim>::
    minimal_temperature (const std::set<types::boundary_id> &) const
    {
      return std::min (inner_temperature, outer_temperature);
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    maximal_temperature (const std::set<types::boundary_id> &) const
    {
      return std::max (inner_temperature, outer_temperature);
    }



    template <int dim>
    void
    SphericalConstant<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary temperature model");
      {
        prm.enter_subsection("Spherical constant");
        {
          prm.declare_entry ("Outer temperature", "0",
                             Patterns::Double (),
                             "Temperature at the outer boundary (lithosphere water/air). Units: K.");
          prm.declare_entry ("Inner temperature", "6000",
                             Patterns::Double (),
                             "Temperature at the inner boundary (core mantle boundary). Units: K.");
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
          inner_temperature = prm.get_double ("Inner temperature");
          outer_temperature = prm.get_double ("Outer temperature");
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
