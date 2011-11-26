//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__boundary_temperature_spherical_constant_h
#define __aspect__boundary_temperature_spherical_constant_h

#include <aspect/boundary_temperature/interface.h>


namespace aspect
{
  namespace BoundaryTemperature
  {
    /**
     * A class that implements a temperature boundary condition for a spherical
     * shell geometry in which the temperature at the inner and outer surfaces
     * (i.e. at the core-mantle and the mantle-lithosphere/atmosphere boundaries)
     * are constant.
     *
     * @ingroup BoundaryTemperatures
     */
    template <int dim>
    class SphericalConstant : public Interface<dim>
    {
      public:
        /**
         * Return the temperature that is to hold on at a particular location on the
         * boundary of the domain. This function returns the constant temperatures
         * read from the parameter file for the inner and outer boundaries.
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

        /**
         * Declare the parameters this class takes through input files.
         * This class declares the inner and outer boundary temperatures.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Temperatures at the inner and outer boundaries.
         */
        double inner_temperature;
        double outer_temperature;
    };
  }
}


#endif
