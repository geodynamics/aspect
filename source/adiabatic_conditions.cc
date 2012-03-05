/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University,
           Wolfgang Bangerth, Texas A&M University,
     Timo Heister, University of Goettingen, 2008-2011 */
/*                                                                */
/*    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


#include <aspect/adiabatic_conditions.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/box.h>

#include <base/std_cxx1x/bind.h>


namespace aspect
{
  template <int dim>
  AdiabaticConditions<dim>::AdiabaticConditions(const GeometryModel::Interface<dim> &geometry_model,
                                                const GravityModel::Interface<dim>  &gravity_model,
                                                const MaterialModel::Interface<dim> &material_model,
                                                const double                         surface_pressure)
    :
    n_points(1000),
    temperatures(n_points, -1),
    pressures(n_points, -1),
    delta_z (geometry_model.maximal_depth() / (n_points-1)),
    geometry_model (geometry_model)
  {
    //TODO: look up real value!
    const double dTdp = 2.5e-8;

    // start with these values: 1200K
    //TODO: use something real
    temperatures[0] = 1200;
    pressures[0] = surface_pressure;

    // now integrate downward using the explicit Euler method for simplicity
    //
    // note: p'(z) = rho(p,T) * |g|
    //       T'(z) = alpha rho |g| T with alpha=1/rho drho/dT
    // TODO: check formulas!
    double z = delta_z;
    for (unsigned int i=1; i<n_points; ++i, z+=delta_z)
      {
        Assert (i < pressures.size(), ExcInternalError());
        Assert (i < temperatures.size(), ExcInternalError());

        const Point<dim> representative_point = geometry_model.representative_point (z);

        // get material parameters and the magnitude of gravity. we assume
        // that gravity always points along the depth direction. this
        // may not strictly be true always but is likely a good enough
        // approximation here.
        const double density = material_model.density(temperatures[i-1], pressures[i-1],
                                                      representative_point);
        const double gravity = gravity_model.gravity_vector(representative_point).norm();

        pressures[i] = pressures[i-1]
                       + density * gravity * delta_z;
        temperatures[i] = temperatures[i-1] +
                          dTdp * density * gravity * delta_z;
      }

    Assert (*min_element (pressures.begin(), pressures.end()) >= 0, ExcInternalError());
    Assert (*min_element (temperatures.begin(), temperatures.end()) >= 0, ExcInternalError());
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
    Assert (i >= 0, ExcInternalError());
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
    Assert (i >= 0, ExcInternalError());
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
  template class AdiabaticConditions<deal_II_dimension>;
}
