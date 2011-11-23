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
#include <aspect/equation_data.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <base/std_cxx1x/bind.h>


namespace aspect
{
  /**
   * A namespace for functions that can convert between an arbitrary point and
   * the depth coordinate. These functions are used in the constructor of the
   * AdiabaticConditions class to store an object that can do this conversion
   * later on without having to remember what exactly the geometry was.
   */
  namespace PointToDepthConversion
  {
    template <int dim>
    double spherical_shell (const Point<dim> &p,
                            const double R0,
                            const double R1)
    {
      // clamp the depth to be positive, can happen due to rounding errors on the mesh.
      // also limit the depth to the maximally possible for this geometry so that
      // we don't get into trouble when looking things up from tables
      return std::min(std::max(R1 - p.norm(),
                               0.0),
                      R1-R0);
    }
  }


  template <int dim>
  AdiabaticConditions<dim>::AdiabaticConditions(const GeometryModel::Interface<dim> &geometry_model,
                                                const GravityModel::Interface<dim>  &gravity_model,
                                                const aspect::MaterialModel::Interface<dim> &material_model)
    :
    n_points(1000),
    temperatures(n_points, -1),
    pressures(n_points, -1)
  {
    // first figure out some global properties of the geometry model
    if (const GeometryModel::SphericalShell<dim> *
        geometry = dynamic_cast<const GeometryModel::SphericalShell<dim>*>(&geometry_model))
      {
        const double R0 = dynamic_cast<const GeometryModel::SphericalShell<dim>&>(geometry_model).inner_radius();
        const double R1 = dynamic_cast<const GeometryModel::SphericalShell<dim>&>(geometry_model).outer_radius();
        delta_z = (R1-R0)/(n_points-1);
        point_to_depth_converter = std_cxx1x::bind (&PointToDepthConversion::spherical_shell<dim>,
                                                    std_cxx1x::_1,
                                                    R0, R1);
      }
    else
      // the following pretty much assumes that we have a spherical shell, for example
      // when we pick a representative point. the code needs to be audited if we
      // have other geometries
      Assert (false, ExcNotImplemented());


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
          = Point<dim>::unit_vector(0) *
            (dynamic_cast<const GeometryModel::SphericalShell<dim>&>(geometry_model).outer_radius()-z);

        const double density = material_model.density(temperatures[i-1], pressures[i-1], representative_point);

        pressures[i] = (pressures[i-1]
                        + pressures[i-1] * 2/z
                        - density *
                        (gravity_model.gravity_vector(representative_point)*Point<dim>::unit_vector(0))
                        * delta_z);
        temperatures[i] = (temperatures[i-1] -
                           dTdp * density *
                           (gravity_model.gravity_vector(representative_point)*Point<dim>::unit_vector(0))
                           * delta_z);
      }

    Assert (*min_element (pressures.begin(), pressures.end()) >= 0, ExcInternalError());
    Assert (*min_element (temperatures.begin(), temperatures.end()) >= 0, ExcInternalError());
  }



  template <int dim>
  double AdiabaticConditions<dim>::pressure (const Point<dim> &p) const
  {
    const double z = point_to_depth_converter(p);

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
    const double z = point_to_depth_converter(p);

    const unsigned int i = z/delta_z;
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
