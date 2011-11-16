//-------------------------------------------------------------
//    $Id: simulator.h 232 2011-10-19 13:30:15Z bangerth $
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__adiabatic_conditions_h
#define __aspect__adiabatic_conditions_h


#include <aspect/material_model_base.h>
#include <deal.II/base/point.h>


namespace aspect
{
  using namespace dealii;


  /**
   * A class that represents adiabatic conditions, i.e. that starts at the top of
   * the domain and integrates pressure and temperature as we go down into depth.
   *
   * @note The implementation has numerous deficiencies indicated in the .cc file
   * and may not quite compute what we want. Specifically, it doesn't currently
   * take into account all the physical parameters it needs, and it also doesn't
   * get gravity right with the exception of the simplest cases.
   */
  template <int dim>
  class AdiabaticConditions
  {
    public:
      /**
       * Constructor. Compute the adiabatic conditions along a vertical
       * transect based on the given material mode.
       */
      AdiabaticConditions (const MaterialModel::Interface<dim> &material_model);

      /**
       * Return the adiabatic temperature at a given point of the domain.
       */
      double temperature (const Point<dim> &p) const;

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
       * depth at which we have computed them. The public member functions
       * of this class interpolate linearly between these points.
       */
      std::vector<double> temperatures, pressures;
  };

}


#endif
