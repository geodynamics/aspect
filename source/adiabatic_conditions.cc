/* $Id: step-32.cc 234 2011-10-19 18:07:35Z bangerth $ */
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
#include <aspect/equation_data.h>


namespace aspect
{
  template <int dim>
  AdiabaticConditions<dim>::AdiabaticConditions(const aspect::MaterialModel::Interface<dim> &material_model)
    :
    n_points(1000),
    temperatures(n_points, -1),
    pressures(n_points, -1)
  {
    //TODO: make sure we use a real geometry model
    const double delta_z = (EquationData::R1-EquationData::R0)/(n_points-1);
    //TODO: look up real value!
    const double dTdp = 2.5e-8;

    // start with these values: 1200K, 1MPa
    //TODO: use something real
    temperatures[0] = 1200;
    pressures[0] = 1e6;

    // now integrate downward using the explicit Euler method for simplicity
    //
    // note: p'(z) = rho(p,T) * g
    //       T'(z) = dT/dp|s dp/dz = dT/dp|S rho(p,T) * g
    double z = delta_z;
    for (unsigned int i=1; i<n_points; ++i, z+=delta_z)
      {
        Assert (i < pressures.size(), ExcInternalError());
        Assert (i < temperatures.size(), ExcInternalError());

        //TODO: use the real gravity model here as a function of z
        const Point<dim> representative_point
          = Point<dim>::unit_vector(0) * (EquationData::R1-z);

        const double density = material_model.density(temperatures[i-1], pressures[i-1], representative_point);

        pressures[i] = (pressures[i-1]
                        + pressures[i-1] * 2/z
                        - density *
                        (EquationData::gravity_vector(representative_point)*Point<dim>::unit_vector(0)) * delta_z);
        temperatures[i] = (temperatures[i-1] -
                           dTdp * density *
                           (EquationData::gravity_vector(representative_point)*Point<dim>::unit_vector(0)) * delta_z);
      }

    Assert (*min_element (pressures.begin(), pressures.end()) >= 0, ExcInternalError());
    Assert (*min_element (temperatures.begin(), temperatures.end()) >= 0, ExcInternalError());
  }



  template <int dim>
  double AdiabaticConditions<dim>::pressure (const Point<dim> &p) const
  {
    const double delta_z = (EquationData::R1-EquationData::R0)/(n_points-1);

    // clamp the depth to be positive, can happen due to rounding errors on the mesh
    const double z = std::min(std::max(EquationData::R1 - p.norm(),
                                       0.0),
                              EquationData::R1-EquationData::R0 - delta_z);

    const unsigned int i = z/delta_z;
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
    const double delta_z = (EquationData::R1-EquationData::R0)/(n_points-1);

    // clamp the depth to be positive, can happen due to rounding errors on the mesh
    const double z = std::min(std::max(EquationData::R1 - p.norm(),
                                       0.0),
                              EquationData::R1-EquationData::R0 - delta_z);

    const unsigned int i = z/delta_z;
    Assert (i >= 0, ExcInternalError());
    Assert (i+1 < temperatures.size(), ExcInternalError());

    // now do the linear interpolation
    const double d=1.0+i-z/delta_z;
    Assert ((d>=0) && (d<=1), ExcInternalError());

    return d*temperatures[i]+(1-d)*temperatures[i+1];
  }
}


//instantiation:
namespace aspect
{
  template class AdiabaticConditions<deal_II_dimension>;
}