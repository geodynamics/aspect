/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/boundary_temperature/spherical_constant.h>
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
      // verify that the geometry is a spherical shell, a chunk, or an
      // ellipsoidal chunk since only for geometries based on spherical shells
      // do we know which boundary indicators are used and what they mean
      const GeometryModel::Interface<dim> *geometry_model = &this->get_geometry_model();
      Assert ((dynamic_cast<const GeometryModel::SphericalShell<dim>*>(geometry_model) != 0
               || dynamic_cast<const GeometryModel::Chunk<dim>*>(geometry_model) != 0
               || dynamic_cast<const GeometryModel::EllipsoidalChunk<dim>*>(geometry_model) != 0),
              ExcMessage ("This boundary model is only implemented if the geometry "
                          "is a spherical shell, ellipsoidal chunk or chunk."));

      const std::string boundary_name = geometry_model->translate_id_to_symbol_name(boundary_indicator);

      if (boundary_name == "inner")
        return inner_temperature;
      else if (boundary_name =="outer")
        return outer_temperature;
      else
        {
          Assert (false, ExcMessage ("Unknown boundary indicator. The given boundary should be ``inner'' or ``outer''."));
          return std::numeric_limits<double>::quiet_NaN();
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
