/*
  Copyright (C) 2011, 2012, 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__geometry_model_interface_h
#define __aspect__geometry_model_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>

#include <set>


namespace aspect
{
  /**
   * A namespace for the definition of properties of the geometry. This
   * primarily includes the definition of the shape of the domain (e.g.
   * whether it is a full spherical shell, a quadrant/octant, a description of
   * the geoid, etc. The classes and functions of this namespace also describe
   * which kinds of boundary conditions hold on the different parts of the
   * boundary of the geometry.
   *
   * @ingroup GeometryModels
   */
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * Base class for classes that describe particular geometries for the
     * domain. These classes must also be able to create coarse meshes and
     * describe what kinds of boundary conditions hold where.
     *
     * @ingroup GeometryModels
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual void initialize ();

        /**
         * Generate a coarse mesh for the geometry described by this class.
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const = 0;

        /**
         * Return the typical length scale one would expect of features in
         * this geometry, assuming realistic parameters.
         *
         * The result of this function is used in computing the scaling factor
         * for the pressure in the Stokes equation. There, the scaling factor
         * is chosen as the ratio of the reference viscosity divided by the
         * length scale. As an example, in the step-32 tutorial program we
         * have determined that a suitable length scale for scaling is 10km,
         * in a domain that is 12,000km across. This length scale suitably
         * matches the order of magnitude for the diameter of plumes in the
         * earth.
         */
        virtual
        double length_scale () const = 0;

        /**
         * Returns the depth that corresponds to the given position. The
         * returned value is between 0 and maximal_depth(), where 0 denotes
         * the surface.
         */
        virtual
        double depth(const Point<dim> &position) const = 0;

        /**
         * Returns a representative point for a given depth. Such a point must
         * lie inside the domain for sure (assuming the given depth is valid),
         * but it is not important where exactly in the domain it is. A good
         * choice would be on a middle axis, for example.
         *
         * The function is used, for example, in computing an initial
         * adiabatic profile. For this, we need to be able to query the
         * density as a function of (adiabatic) pressure and temperature.
         * However, the functions returning the density from the material
         * model also depend on location. Since we are interested only in
         * computing a depth-dependent adiabatic profile without lateral
         * variation, we need to be able to query the density at
         * "representative points" at given depths, without caring too much
         * where exactly that is -- but at points that we know for sure are
         * inside the domain.
         */
        virtual
        Point<dim> representative_point(const double depth) const = 0;

        /**
         * Returns the maximal depth of this geometry.
         */
        virtual
        double maximal_depth() const = 0;


        /**
         * Return the set of boundary indicators that are used by this model.
         * This information is used to determine what boundary indicators can
         * be used in the input file.
         */
        virtual
        std::set<types::boundary_id>
        get_used_boundary_indicators () const = 0;

        /**
         * Return a mapping from symbolic names of each part of the boundary
         * to the corresponding boundary indicator. This allows users to
         * specify *names*, not just *numbers* in their input files when
         * describing which parts of the boundary have to satisfy which
         * boundary conditions.
         *
         * An example would be that the "box" geometry returns a map of the
         * form <code>{{"left"->0}, {"right"->1}, {"bottom"->2},
         * {"top"->3}}</code> in 2d.
         *
         * The default implementation of this function returns an empty map.
         * This still allows the use of a geometry model that does not
         * implement this function but forces the user to identify parts of
         * the boundary by their boundary indicator number, rather than using
         * a symbolic name.
         *
         * @note Names may contain spaces and numbers, but they may not
         * contain special characters and they should not equal the text
         * representation of numbers (e.g., a name "10" is ill-advised).
         *
         * @note Since in practice boundary indicators can be provided either
         * via number or symbolic name, the mapping from something given in
         * the input is not entirely trivial -- in particular, because a
         * function also has to do some error checking that a given string in
         * fact matches any known boundary indicator. To this end, use
         * GeometryModel::translate_boundary_indicator() and
         * GeometryModel::translate_boundary_indicators().
         *
         * @return A map from symbolic names to boundary indicators. The map
         * should provide a symbolic name for each used boundary indicator as
         * returned by get_user_boundary_indicators(). In the end, however, a
         * geometry model may define multiple symbolic names for the same
         * boundary or not define any.
         */
        virtual
        std::map<std::string,types::boundary_id>
        get_symbolic_boundary_names_map () const;

        /**
         * For a given name of a boundary component, translate it to its
         * numeric value -- either by using one of the symbolic values in the
         * mapping that derived classes provide through the
         * get_symbolic_boundary_names_map() function, or by converting its
         * string representation into a number.
         *
         * @param name A name or number (as string). Leading and trailing
         * spaces are removed from this string.
         * @return A boundary indicator number corresponding to the given
         * name. If the name does not represent either a symbolic name or a
         * number, this function will throw an exception of type std::string
         * that explains the error.
         */
        types::boundary_id
        translate_symbolic_boundary_name_to_id (const std::string &name) const;

        /**
         * For each one of the given names of boundary components, translate
         * it to its numeric value -- either by using one of the symbolic
         * values in the mapping that derived classes provide through the
         * get_symbolic_boundary_names_map() function, or by converting its
         * string representation into a number.
         *
         * @param names A list of names or numbers (as strings). Leading and
         * trailing spaces are removed from the strings.
         * @return A list of boundary indicator numbers corresponding to the
         * given list of names. If one of the given names does not represent
         * either a symbolic name or a number, this function will throw an
         * exception of type std::string that explains the error.
         */
        std::vector<types::boundary_id>
        translate_symbolic_boundary_names_to_ids (const std::vector<std::string> &names) const;

        /**
         * Given a boundary indicator, try and see whether this geometry model
         * provides a symbolic name for it. If the geometry model does not
         * provide a symbolic name, return the empty string.
         *
         * @param boundary_id The boundary indicator to be translated.
         * @return A string representation for this boundary indicator, if one
         * is available.
         */
        std::string
        translate_id_to_symbol_name (const types::boundary_id boundary_id) const;

        /**
         * Returns a set of periodic boundary pairs.  The elements of the set
         * are a pair of boundary ids and a cartesian unit direction each. The
         * base class returns an empty set, so this does nothing unless you
         * specifically use a geometry model with periodic boundary conditions
         */
        virtual
        std::set< std::pair< std::pair<types::boundary_id, types::boundary_id>, unsigned int> >
        get_periodic_boundary_pairs () const;

        /**
         * If true, the geometry contains cells with boundaries that are not straight
         * and have a deal.II boundary object attached to it. If the return value is
         * @p false, certain operation can be optimized.The default implementation of
         * this function will return @p true.
         */
        virtual
        bool
        has_curved_elements() const;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

    };



    /**
     * Register a geometry model so that it can be selected from the parameter
     * file.
     *
     * @param name A string that identifies the geometry model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this geometry model wants to read
     * from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this geometry model.
     *
     * @ingroup GeometryModels
     */
    template <int dim>
    void
    register_geometry_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> *(*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The geometry model will also be asked to read its runtime parameters
     * already.
     *
     * @ingroup GeometryModels
     */
    template <int dim>
    Interface<dim> *
    create_geometry_model (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered geometry models.
     *
     * @ingroup GeometryModels
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a geometry model, register it with the functions that can declare
     * their parameters and create these objects.
     *
     * @ingroup GeometryModels
     */
#define ASPECT_REGISTER_GEOMETRY_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_GEOMETRY_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::GeometryModel::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::GeometryModel::register_geometry_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::GeometryModel::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::GeometryModel::register_geometry_model<3>, \
                                name, description); \
  }
  }
}


#endif
