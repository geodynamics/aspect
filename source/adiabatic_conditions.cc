/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/adiabatic_conditions.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/box.h>

#include <deal.II/base/std_cxx1x/bind.h>


namespace aspect
{
  template <int dim>
  AdiabaticConditions<dim>::AdiabaticConditions(const GeometryModel::Interface<dim> &geometry_model,
                                                const GravityModel::Interface<dim>  &gravity_model,
                                                const MaterialModel::Interface<dim> &material_model,
                                                const CompositionalInitialConditions::Interface<dim> &compositional_initial_conditions,
                                                const double                         surface_pressure,
                                                const double                         surface_temperature,
                                                const unsigned int                   n_compositional_fields)
    :
    n_points(1000),
    temperatures(n_points, -1),
    pressures(n_points, -1),
    delta_z (geometry_model.maximal_depth() / (n_points-1)),
    geometry_model (geometry_model)
  {

    temperatures[0] = surface_temperature;
    pressures[0]    = surface_pressure;

    // now integrate downward using the explicit Euler method for simplicity
    //
    // note: p'(z) = rho(p,T) * |g|
    //       T'(z) = alpha |g| T / c_p
    double z;
    for (unsigned int i=1; i<n_points; ++i)
      {
        Assert (i < pressures.size(), ExcInternalError());
        Assert (i < temperatures.size(), ExcInternalError());

        z = double(i)/double(n_points-1)*geometry_model.maximal_depth();

        const Point<dim> representative_point = geometry_model.representative_point (z);

        typename MaterialModel::Interface<dim>::MaterialModelInputs in(1, n_compositional_fields);
        typename MaterialModel::Interface<dim>::MaterialModelOutputs out(1, n_compositional_fields);
        in.position[0] = representative_point;
        in.temperature[0] = temperatures[i-1];
        in.pressure[0] = pressures[i-1];

        //TODO: we look up the composition at the representative point, but we should
        // use averaged compositional values here. Right?
        for (unsigned int c=0; c<n_compositional_fields; ++c)
          in.composition[0][c] = compositional_initial_conditions.initial_composition(representative_point, c);

        in.strain_rate.resize(0); // we do not need the viscosity
        material_model.evaluate(in, out);

        // get the magnitude of gravity. we assume
        // that gravity always points along the depth direction. this
        // may not strictly be true always but is likely a good enough
        // approximation here.
        const double density = out.densities[0];
        const double alpha = out.thermal_expansion_coefficients[0];
        const double cp = out.specific_heat[0];
        const double gravity = gravity_model.gravity_vector(representative_point).norm();

        pressures[i] = pressures[i-1]
                       + density * gravity * delta_z;
        temperatures[i] = temperatures[i-1] * (1 +
                                               alpha * gravity * delta_z / cp);
      }

    Assert (*std::min_element (pressures.begin(), pressures.end()) >=
            -std::numeric_limits<double>::epsilon() * pressures.size(),
            ExcInternalError());
    Assert (*std::min_element (temperatures.begin(), temperatures.end()) >=
            -std::numeric_limits<double>::epsilon() * temperatures.size(),
            ExcInternalError());
  }

  template <int dim>
  void AdiabaticConditions<dim>::get_adiabatic_temperature_profile(std::vector<double> &values) const
  {
    const unsigned int num_slices = values.size();
    const double max_depth = geometry_model.maximal_depth();
    double depth = 0.0;

    for (unsigned int n = 0 ; n < num_slices; n++)
      {
        depth = n * max_depth / (num_slices-1);
        const Point<dim> p = geometry_model.representative_point(depth);
        values[n] = temperature(p);
      }
  }


  template <int dim>
  double AdiabaticConditions<dim>::pressure (const Point<dim> &p) const
  {
    const double z = geometry_model.depth(p);

    if (z >= geometry_model.maximal_depth())
      {
        Assert (z <= geometry_model.maximal_depth() + delta_z,
                ExcInternalError());
        return pressures.back();
      }

    const unsigned int i = static_cast<unsigned int>(z/delta_z);
    Assert ((z/delta_z) >= 0, ExcInternalError());
    Assert (i+1 < pressures.size(), ExcInternalError());

    // now do the linear interpolation
    const double d=1.0+i-z/delta_z;
    Assert ((d>=0) && (d<=1), ExcInternalError());

    return d*pressures[i]+(1-d)*pressures[i+1];
  }



  template <int dim>
  double AdiabaticConditions<dim>::temperature (const Point<dim> &p) const
  {
    const double z = geometry_model.depth(p);

    if (z >= geometry_model.maximal_depth())
      {
        Assert (z <= geometry_model.maximal_depth() + delta_z,
                ExcInternalError());
        return temperatures.back();
      }

    const unsigned int i = static_cast<unsigned int>(z/delta_z);
    Assert ((z/delta_z) >= 0, ExcInternalError());
    Assert (i+1 < temperatures.size(), ExcInternalError());

    // now do the linear interpolation
    const double d=1.0+i-z/delta_z;
    Assert ((d>=0) && (d<=1), ExcInternalError());

    return d*temperatures[i]+(1-d)*temperatures[i+1];
  }
}


// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template class AdiabaticConditions<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
