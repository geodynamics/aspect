//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__geometry_model_spherical_shell_h
#define __aspect__geometry_model_spherical_shell_h

#include <aspect/geometry_model/interface.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class that describes the geometry as a spherical shell. To be more
     * precise, at least in 2d this class can also simulate just a sector
     * of the spherical shell geometry, in particular a half ring and a
     * quarter ring.
     *
     * The parameters that describe this geometry and that are read from
     * the input file are the inner and outer radii of the shell and the
     * opening angle of the section of the shell we want to build.
     */
    template <int dim>
    class SphericalShell : public Interface<dim>
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
         * As discussed in the step-32 tutorial program, an appropriate length scale for
         * this geometry is 10km, so we return $10^4$ here.
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
         * corresponding to the inner and outer parts of the boundary.
         */
        virtual
        std::set<unsigned char>
        get_temperature_dirichlet_boundary_indicators () const;

        /**
         * Declare the parameters this class takes through input files.
         * The default implementation of this function does not describe
         * any parameters. Consequently, derived classes do not have to
         * overload this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file. The default implementation of this function does not read
         * any parameters. Consequently, derived classes do not have to
         * overload this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Return the inner radius of the shell.
         */
        double
        inner_radius () const;

        /**
         * Return the outer radius of the shell.
         */
        double
        outer_radius () const;

        /**
         * Return the opening angle of the shell sector.
         */
        double
        opening_angle () const;

      public:
        /**
         * Inner and outer radii of the spherical shell.
         */
        double R0, R1;

        /**
         * Opening angle of the section of the shell that we simulate.
         */
        double phi;
    };
  }
}


#endif
