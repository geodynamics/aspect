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
/*  $Id$  */


#ifndef __aspect__adiabatic_conditions_h
#define __aspect__adiabatic_conditions_h


#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/compositional_initial_conditions/interface.h>
#include <deal.II/base/point.h>


namespace aspect
{
  using namespace dealii;


  /**
   * A class that represents adiabatic conditions, i.e. that starts at the top
   * of the domain and integrates pressure and temperature as we go down into
   * depth.
   *
   * @note The implementation has numerous deficiencies indicated in the .cc
   * file and may not quite compute what we want. Specifically, it doesn't
   * currently take into account all the physical parameters it needs, and it
   * also doesn't get gravity right with the exception of the simplest cases.
   */
  template <int dim>
  class AdiabaticConditions
  {
    public:
      /**
       * Constructor. Compute the adiabatic conditions along a vertical
       * transect of the geometry based on the given material model and other
       * quantities.
       */
      AdiabaticConditions (const GeometryModel::Interface<dim> &geometry_model,
                           const GravityModel::Interface<dim>  &gravity_model,
                           const MaterialModel::Interface<dim> &material_model,
                           const CompositionalInitialConditions::Interface<dim> &compositional_initial_conditions,
                           const double                         surface_pressure,
                           const double                         surface_temperature,
                           const unsigned int                   n_compositional_fields);

      /**
       * Return the adiabatic temperature at a given point of the domain.
       */
      double temperature (const Point<dim> &p) const;

      /**
       * Return the adiabatic temperature profile as a vector of values
       * corresponding to increasing depth.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void get_adiabatic_temperature_profile(std::vector<double> &values) const;

      /**
       * Return the adiabatic pressure at a given point of the domain.
       */
      double pressure (const Point<dim> &p) const;

    private:
      /**
       * Number of points at which we compute the adiabatic values.
       */
      const unsigned int n_points;

      /**
       * Vectors of values of temperatures and pressures on a transect into
       * depth at which we have computed them. The public member functions of
       * this class interpolate linearly between these points.
       */
      std::vector<double> temperatures, pressures;

      /**
       * Interval spacing between each two data points in the tables above
       * with regard to the depth coordinate.
       */
      double delta_z;

      /**
       * A reference to the geometry model which we need when converting
       * between arbitrary points at which temperature and pressure are
       * interpolated and the depth coordinate we use to pre-compute values.
       */
      const GeometryModel::Interface<dim> &geometry_model;
  };

}


#endif
