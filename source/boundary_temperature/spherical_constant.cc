// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <deal.II/base/signaling_nan.h>
#include <aspect/boundary_temperature/spherical_constant.h>
#include <aspect/geometry_model/sphere.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
    template <int dim>
    void
    SphericalConstant<dim>::
    initialize ()
    {
      // verify that the geometry is supported by this plugin
      AssertThrow ( Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>>(this->get_geometry_model()) ||
                    Plugins::plugin_type_matches<const GeometryModel::EllipsoidalChunk<dim>>(this->get_geometry_model()),
                    ExcMessage ("This boundary model is only implemented if the geometry is "
                                "one of the spherical geometries."));

      // no inner boundary in a full sphere
      if (Plugins::plugin_type_matches<const GeometryModel::Sphere<dim>>(this->get_geometry_model()))
        inner_boundary_indicator = numbers::invalid_unsigned_int;
      else
        inner_boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id("bottom");

      outer_boundary_indicator = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");
    }



    template <int dim>
    double
    SphericalConstant<dim>::
    boundary_temperature (const types::boundary_id boundary_indicator,
                          const Point<dim> &) const
    {
      if (boundary_indicator == outer_boundary_indicator)
        return outer_temperature;
      else if (boundary_indicator == inner_boundary_indicator)
        return inner_temperature;
      else
        AssertThrow (false,
                     ExcMessage ("Unknown boundary indicator for geometry model. "
                                 "The given boundary should be ``top'' or ``bottom''."));

      return numbers::signaling_nan<double>();
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
          prm.declare_entry ("Outer temperature", "0.",
                             Patterns::Double (),
                             "Temperature at the outer boundary (lithosphere water/air). Units: $\\text{K}$.");
          prm.declare_entry ("Inner temperature", "6000.",
                             Patterns::Double (),
                             "Temperature at the inner boundary (core mantle boundary). Units: $\\text{K}$.");
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
