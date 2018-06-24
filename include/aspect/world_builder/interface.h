/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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


#ifndef _aspect_world_builder_interface_h
#define _aspect_world_builder_interface_h

#include <aspect/gravity_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/plugins.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <boost/algorithm/string.hpp>
#include <boost/property_tree/ptree.hpp>

#include <set>


namespace aspect
{
  //template <int dim> class SimulatorAccess;
  /**
   * @ingroup WorldBuilder
   */
  namespace WorldBuilder
  {
    using namespace dealii;

    class WorldBuilderParameterHandler
    {
      public:
        WorldBuilderParameterHandler(const char path_separator = '.');
        void set_tree(boost::property_tree::ptree &tree);
        boost::property_tree::ptree &get_tree();
        void enter_subsection (const std::string &subsection);
        void leave_subsection ();
        std::string get_current_path () const;
        std::string get_current_full_path(const std::string &name) const;

        std::string get(const std::string &entry_string) const;
        std::vector<std::string> get_array(const std::string &s);
        std::vector<std::vector<std::string> > get_double_array(const std::string &s);

        void declare_entry (const std::string           &entry,
                            const std::string           &default_value,
                            const Patterns::PatternBase &pattern,
                            const std::string           &documentation);

        /**
         * A list of patterns that are used to describe the parameters of this
         * object. Every nodes in the property tree corresponding to a parameter
         * stores an index into this array.
         */
        std::vector<std::unique_ptr<const Patterns::PatternBase> > patterns;

      private:
        boost::property_tree::ptree wb_tree;
        std::vector<std::string> subsection_path;
        const char path_separator;

        /**
         * Exception
         */
        DeclException2 (ExcValueDoesNotMatchPattern,
                        std::string, std::string,
                        << "The string <" << arg1
                        << "> does not match the given pattern <" << arg2 << ">.");
    };


    using namespace dealii;

    /**
     * Base class for classes that describe particular initial conditions.
     * All plugins, which describe a plate tectonic features in terms of
     * temperature, composition etc, are derived from this class.
     *
     * @ingroup WorldBuilder
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
         * Return the temperature of an object.
         */
        virtual
        double temperature (const Utilities::NaturalCoordinate<3> &position,
                            const unsigned int object_number,
                            const double depth,
                            double temperature) const = 0;

        /**
         * Return the composition of an object.
         */
        virtual
        double composition (const Utilities::NaturalCoordinate<3> &position,
                            const unsigned int n_comp,
                            const unsigned int object_number,
                            const double depth,
                            double composition) const = 0;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (WorldBuilderParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);



      private:


    };



    /**
     * Manager class
     */
    template <int dim>
    class Manager : public SimulatorAccess<dim>,
      public InitialTemperature::Interface<dim>,
      public InitialComposition::Interface<dim>
    {
      public:
        /**
         * The constructor of the world builder manager.
         */
        Manager();

        /**
         * This initializes all the modules/plugins.
         */
        void initialize();


        /**
         * Register a world builder module so that it can be selected from the parameter
         * file.
         *
         * @param name A string that identifies the world builder module
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this world builder module wants
         * to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this world builder module.
         *
         * @ingroup WorldBuilder
         */
        static
        void register_world_builder (const std::string &name,
                                     const std::string &description,
                                     void (*declare_parameters_function) (ParameterHandler &),
                                     Interface<dim> *(*factory_function) ());


        virtual
        double initial_temperature (const dealii::Point<dim> &position) const;

        virtual
        double initial_composition (const dealii::Point<dim> &position, const unsigned int n_comp) const;

        /**
         * Declare the runtime parameters of the registered world builder
         * models.
         *
         * @ingroup WorldBuilder
         */
        void
        declare_parameters (WorldBuilderParameterHandler &prm);

        /**
         * Parse the runtime parameters of the registered world builder
         * models.
         *
         * @ingroup WorldBuilder
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Declare the runtime parameters of the registered world builder
         * models.
         *
         * @ingroup WorldBuilder
         */
        static void
        declare_parameters (ParameterHandler &prm);

        ParameterHandler wb_prm;
        static WorldBuilderParameterHandler ph;
        std::vector<std_cxx11::shared_ptr<Interface<dim> > >  modules;
        std::vector<std_cxx11::shared_ptr<Interface<dim> > >  ordered_modules;
        std::vector<std::string> module_names;

      private:
        double rotation_angle;
        bool has_sticky_air;
        double sticky_air_thickness;
        unsigned int sticky_air_composition;
        bool has_top_crust;
        double top_crust_tickness;
        unsigned int top_crust_composition;
        std::vector<double> rotation_point;
        double surface_temperature;
        unsigned int number_of_objects;

        double potential_mantle_temperature, thermal_expansion_coefficient_alfa, specific_heat_Cp;

        //const std_cxx11::shared_ptr<GravityModel::Interface<dim> > gravity_model;
        //const std_cxx11::shared_ptr<GeometryModel::Interface<dim> > geometry_model;

    };



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a initial topography model, register it with the functions that
     * can declare their parameters and create these objects.
     *
     * @ingroup WorldBuilder
     */
#define ASPECT_REGISTER_WORLD_BUILDER(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_WORLD_BUILDER ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::WorldBuilder::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::WorldBuilder::Manager<2>::register_world_builder, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::WorldBuilder::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::WorldBuilder::Manager<3>::register_world_builder, \
                                name, description); \
  }

#define ASPECT_REGISTER_WORLD_BUILDER_AS_INITIAL_TEMPERATURE_MODEL(classname,name,description) \
  namespace ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialTemperature::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::InitialTemperature::Manager<2>::register_initial_temperature, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialTemperature::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::InitialTemperature::Manager<3>::register_initial_temperature, \
                                name, description); \
  }

#define ASPECT_REGISTER_WORLD_BUILDER_AS_INITIAL_COMPOSITION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialComposition::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::InitialComposition::Manager<2>::register_initial_composition, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::InitialComposition::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::InitialComposition::Manager<3>::register_initial_composition, \
                                name, description); \
  }
  }

}


#endif
