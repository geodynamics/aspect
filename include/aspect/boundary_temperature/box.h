//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__boundary_temperature_box_h
#define __aspect__boundary_temperature_box_h

#include <aspect/boundary_temperature/interface.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
    /**
     * A class that implements a temperature boundary condition for a box
     * geometry.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class Box : public Interface<dim>
    {
      public:
        /**
         * Return the temperature that is to hold at a particular location on the
         * boundary of the domain. This function returns constant temperatures
         * at the left and right boundaries.
         *
         * @param geometry_model The geometry model that describes the domain. This may
         *   be used to determine whether the boundary temperature model is implemented
         *   for this geometry.
         * @param boundary_indicator The boundary indicator of the part of the boundary
         *   of the domain on which the point is located at which we are requesting the
         *   temperature.
         * @param location The location of the point at which we ask for the temperature.
         **/
        virtual
        double temperature (const GeometryModel::Interface<dim> &geometry_model,
                            const unsigned int                   boundary_indicator,
                            const Point<dim>                    &location) const;

        /**
         * Return the minimal the temperature on that part of the boundary
         * on which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double minimal_temperature () const;

        /**
         * Return the maximal the temperature on that part of the boundary
         * on which Dirichlet conditions are posed.
         *
         * This value is used in computing dimensionless numbers such as the
         * Nusselt number indicating heat flux.
         */
        virtual
        double maximal_temperature () const;
    };
  }
}


#endif
