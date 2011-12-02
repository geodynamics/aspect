//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__geometry_model_box_h
#define __aspect__geometry_model_box_h

#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes a box geometry.
     */
    template <int dim>
    class Box : public Interface<dim>
    {
      public:
        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

        /**
         * Return the typical length scale one would expect of features in this geometry,
         * assuming realistic parameters.
         *
         * We return 1/100th of the diameter of the box.
         */
        virtual
        double length_scale () const;

        /**
         * Return a set of boundary indicators that correspond to parts of the
         * boundary on which the temperature is fixed, i.e., on which Dirichlet
         * boundary conditions are posed. These boundary conditions for the temperature
         * are then described by classes derived from BoundaryTemperature::Interface.
         *
         * For the geometry used by this class, the returned set is $\{0,1\}$,
         * corresponding to the left and right boundaries.
         */
        virtual
        std::set<unsigned char>
        get_temperature_dirichlet_boundary_indicators () const;

        /**
         * Return a set of boundary indicators that correspond to parts of the
         * boundary on which the velocity is zero.
         *
         * For the geometry here, we prescribe zero velocity on all
         * boundaries on which the temperature is not prescribed.
         */
        virtual
        std::set<unsigned char>
        get_zero_velocity_boundary_indicators () const;

        /**
         * Return a set of boundary indicators that correspond to parts of the
         * boundary on which the velocity is tangential but may not be zero.
         *
         * For the present geometry, we prescribe tangential flow on
        * all boundaries where the temperature is prescribed by the
        * get_temperature_dirichlet_boundary_indicators() function.
         */
        virtual
        std::set<unsigned char>
        get_tangential_velocity_boundary_indicators () const;
    };
  }
}


#endif
