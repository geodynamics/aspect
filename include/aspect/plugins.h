/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_plugins_h
#define _aspect_plugins_h

#include <aspect/global.h>

#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <tuple>
#include <deal.II/base/exceptions.h>

#include <boost/core/demangle.hpp>

#include <string>
#include <list>
#include <set>
#include <map>
#include <iostream>
#include <typeinfo>


namespace aspect
{
  template <int dim> class SimulatorAccess;

  namespace Plugins
  {
    using namespace dealii;

    /**
     * This function returns if a given plugin (e.g. a material model returned
     * from SimulatorAccess::get_material_model() ) matches a certain plugin
     * type (e.g. MaterialModel::Simple). This check is needed, because often
     * it is only possible to get a reference to an Interface, not the actual
     * plugin type, but the actual plugin type might be important. For example
     * a radial gravity model might only be implemented for spherical geometry
     * models, and would want to check if the current geometry is in fact a
     * spherical shell.
     */
    template <typename TestType, typename PluginType>
    inline
    bool
    plugin_type_matches (const PluginType &object)
    {
      return (dynamic_cast<const TestType *> (&object) != nullptr);
    }

    /**
     * This function converts a reference to a type (in particular a reference
     * to an interface class) into a reference to a different type (in
     * particular a plugin class). This allows accessing members of the plugin
     * that are not specified in the interface class. Note that you should
     * first check if the plugin type is actually convertible by calling
     * plugin_matches_type() before calling this function. If the plugin is
     * not convertible this function throws an exception.
     */
    template <typename TestType, typename PluginType>
    inline
    TestType &
    get_plugin_as_type (PluginType &object)
    {
      AssertThrow(plugin_type_matches<TestType>(object),
                  ExcMessage("You have requested to convert a plugin of type <"
                             + boost::core::demangle(typeid(PluginType).name())
                             + "> into type <"
                             + boost::core::demangle(typeid(TestType).name()) +
                             ">, but this cast cannot be performed."));

      // We can safely dereference the pointer, because we checked above that
      // the object is actually of type TestType, and so the result
      // is not a nullptr.
      return *dynamic_cast<TestType *> (&object);
    }
  }

  namespace internal
  {
    /**
     * A namespace for the definition of classes that have to do with the
     * plugin architecture of ASPECT.
     */
    namespace Plugins
    {
      using namespace dealii;

      /**
       * An internal class that is used in the definition of the
       * ASPECT_REGISTER_* macros. Given a registration function, a classname,
       * a description of what it does, and a name for the parameter file, it
       * registers the model with the proper authorities.
       *
       * The registration happens in the constructor. The typical use case of
       * this function is thus the creation of a dummy object in some
       * otherwise unused namespace.
       */
      template <typename InterfaceClass,
                typename ModelClass>
      struct RegisterHelper
      {
        /**
         * Constructor. Given a pointer to a registration function and name
         * and description of the class, this constructor registers the class
         * passed as second template argument.
         */
        RegisterHelper (void (*register_function) (const std::string &,
                                                   const std::string &,
                                                   void ( *)(ParameterHandler &),
                                                   InterfaceClass * ( *)()),
                        const char *name,
                        const char *description)
        {
          register_function (name,
                             description,
                             &ModelClass::declare_parameters,
                             &factory);
        }

        /**
         * A factory object that just creates object of the type registered by
         * this class.
         */
        static
        InterfaceClass *factory ()
        {
          return new ModelClass();
        }
      };


      /**
       * A class that stores a list of registered plugins for the given
       * interface type.
       */
      template <typename InterfaceClass>
      struct PluginList
      {
        /**
         * A type describing everything we need to know about a plugin.
         *
         * The entries in the tuple are: - The name by which it can be
         * selected. - A description of this plugin that will show up in the
         * documentation in the parameter file. - A function that can declare
         * the run-time parameters this plugin takes from the parameter file.
         * - A function that can produce objects of this plugin type.
         */
        using PluginInfo
        = std::tuple<std::string,
        std::string,
        void ( *) (ParameterHandler &),
        InterfaceClass *( *) ()>;

        /**
         * A pointer to a list of all registered plugins.
         *
         * The object is a pointer rather than an object for the following
         * reason: objects with static initializers (such as =0) are
         * initialized before any objects for which one needs to run
         * constructors. consequently, we can be sure that this pointer is set
         * to zero before we ever try to register a postprocessor, and
         * consequently whenever we run Manager::register_postprocessor, we
         * need not worry whether we try to add something to this list before
         * the lists's constructor has successfully run
         */
        static std::list<PluginInfo> *plugins;

        /**
         * Destructor.
         */
        ~PluginList ();

        /**
         * Register a plugin by name, description, parameter declaration
         * function, and factory function. See the discussion for the
         * PluginInfo type above for more information on their meaning.
         */
        static
        void register_plugin (const std::string &name,
                              const std::string &description,
                              void (*declare_parameters_function) (ParameterHandler &),
                              InterfaceClass * (*factory_function) ());

        /**
         * Generate a list of names of the registered plugins separated by '|'
         * so that they can be taken as the input for Patterns::Selection.
         *
         * To make it easier to visually scan through the list of plugins,
         * names are sorted alphabetically.
         */
        static
        std::string get_pattern_of_names ();

        /**
         * Return a string that describes all registered plugins using the
         * descriptions that have been provided at the time of registration.
         *
         * To make it easier to visually scan through the list of plugins,
         * names are sorted alphabetically.
         */
        static
        std::string get_description_string ();

        /**
         * Let all registered plugins define their parameters.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Given the name of one plugin, create a corresponding object and
         * return a pointer to it. The second argument provides a hint where
         * this function was called from, to be printed in case there is an
         * error.
         *
         * Ownership of the object is handed over to the caller of this
         * function.
         */
        static
        InterfaceClass *
        create_plugin (const std::string  &name,
                       const std::string &documentation);

        /**
         * Given the name of one plugin, create a corresponding object and
         * return a pointer to it. The second argument provides a hint where
         * this function was called from, to be printed in case there is an
         * error. Before returning, let the newly created object read its run-
         * time parameters from the parameter object.
         *
         * Ownership of the object is handed over to the caller of this
         * function.
         */
        static
        InterfaceClass *
        create_plugin (const std::string  &name,
                       const std::string &documentation,
                       ParameterHandler &prm);

        /**
         * For the current plugin subsystem, write a connection graph of all of the
         * plugins we know about, in the format that the
         * programs dot and neato understand. This allows for a visualization of
         * how all of the plugins that ASPECT knows about are interconnected, and
         * connect to other parts of the ASPECT code.
         *
         * @param plugin_system_name The name to be used for the current
         *   plugin system. This name will be used for the "Interface"
         *   class to which all plugins connect.
         * @param output_stream The stream to write the output to.
         * @param attachment_point The point to which a plugin subsystem
         *   feeds information. By default, this is the Simulator class,
         *   but some plugin systems (most notably the visualization
         *   postprocessors, which feeds to one of the postprocessor
         *   classes) hook into other places. If other than the
         *   "Simulator" default, the attachment point should be of
         *   the form <code>typeid(ClassName).name()</code> as this is
         *   the form used by this function to identify nodes in the
         *   plugin graph.
         */
        static
        void
        write_plugin_graph (const std::string &plugin_system_name,
                            std::ostream      &output_stream,
                            const std::string &attachment_point = "Simulator");

        /**
         * Exception.
         */
        DeclException1 (ExcUnknownPlugin,
                        std::string,
                        << "Can't create a plugin of name <" << arg1
                        << "> because such a plugin hasn't been declared.");
      };


      /* ------------------------ template and inline functions --------------------- */

      template <typename InterfaceClass>
      PluginList<InterfaceClass>::
      ~PluginList ()
      {
        // if any plugins have been registered, then delete
        // the list
        if (plugins != nullptr)
          delete plugins;
        plugins = nullptr;
      }



      template <typename InterfaceClass>
      void
      PluginList<InterfaceClass>::
      register_plugin (const std::string &name,
                       const std::string &description,
                       void (*declare_parameters_function) (ParameterHandler &),
                       InterfaceClass * (*factory_function) ())
      {
        // see if this is the first time we get into this
        // function and if so initialize the static member variable
        if (plugins == nullptr)
          plugins = new std::list<PluginInfo>();

        // verify that the same name has not previously been
        // used to register a plugin, since we would then no
        // longer be able to identify the plugin
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          Assert (std::get<0>(*p) != name,
                  ExcMessage ("A plugin with name <" + name + "> has "
                              "already been registered!"));

        // now add one record to the list
        plugins->push_back (PluginInfo(name,
                                       description,
                                       declare_parameters_function,
                                       factory_function));
      }



      template <typename InterfaceClass>
      std::string
      PluginList<InterfaceClass>::
      get_pattern_of_names ()
      {
        Assert (plugins != nullptr,
                ExcMessage ("No plugins registered!?"));

        // get all names and put them into a data structure that keeps
        // them sorted
        std::set<std::string> names;
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          names.insert (std::get<0>(*p));

        // now create a pattern from all of these sorted names
        std::string pattern_of_names;
        for (typename std::set<std::string>::const_iterator
             p = names.begin();
             p != names.end(); ++p)
          {
            if (pattern_of_names.size() > 0)
              pattern_of_names += "|";
            pattern_of_names += *p;
          }

        return pattern_of_names;
      }



      template <typename InterfaceClass>
      std::string
      PluginList<InterfaceClass>::
      get_description_string ()
      {
        std::string description;

        // get all names_and_descriptions and put them into a data structure that keeps
        // them sorted
        std::map<std::string,std::string> names_and_descriptions;
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          names_and_descriptions[std::get<0>(*p)] = std::get<1>(*p);;

        // then output it all
        std::map<std::string,std::string>::const_iterator
        p = names_and_descriptions.begin();
        while (true)
          {
            // write the name and
            // description of the
            // parameter
            description += "`";
            description += p->first;
            description += "': ";
            description += p->second;

            // increment the pointer
            // by one. if we are not
            // at the end yet then
            // add an empty line
            ++p;
            if (p != names_and_descriptions.end())
              description += "\n\n";
            else
              break;
          }

        return description;
      }




      template <typename InterfaceClass>
      void
      PluginList<InterfaceClass>::
      declare_parameters (ParameterHandler &prm)
      {
        Assert (plugins != nullptr,
                ExcMessage ("No postprocessors registered!?"));

        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          (std::get<2>(*p))(prm);
      }



      template <typename InterfaceClass>
      InterfaceClass *
      PluginList<InterfaceClass>::
      create_plugin (const std::string &name,
                     const std::string &documentation)
      {
        (void)documentation;
        Assert (plugins != nullptr,
                ExcMessage ("No postprocessors registered!?"));
        AssertThrow (name != "unspecified",
                     ExcMessage(std::string("A plugin must have a name!\n\n"
                                            "This function was asked to create a plugin but no name for the "
                                            "plugin was provided. This may be due to the fact that you did not "
                                            "explicitly specify a name for this plugin in your input file and "
                                            "ASPECT does not provide a default for this kind of plugin, for "
                                            "example because no generally useful plugin exists. An example "
                                            "is that there is no default geometry: You need to explicitly "
                                            "provide one in the input file, and it seems like you have not "
                                            "done so.\n\n"
                                            "To find out which kind of plugin this function tries to create, "
                                            "take a look at the backtrace of this error message.\n\n"
                                            "The place that called this function also provided as "
                                            "additional information this:\n\n"
                                            "   <")
                                + documentation + ">"));

        for (typename std::list<PluginInfo>::const_iterator p = plugins->begin();
             p != plugins->end(); ++p)
          if (std::get<0>(*p) == name)
            {
              InterfaceClass *i = std::get<3>(*p)();
              return i;
            }

        AssertThrow (false, ExcUnknownPlugin(name));
        return nullptr;
      }



      template <typename InterfaceClass>
      InterfaceClass *
      PluginList<InterfaceClass>::
      create_plugin (const std::string &name,
                     const std::string &documentation,
                     ParameterHandler  &prm)
      {
        InterfaceClass *i = create_plugin(name, documentation);
        i->parse_parameters (prm);
        return i;
      }



      template <typename InterfaceClass>
      void
      PluginList<InterfaceClass>::
      write_plugin_graph (const std::string &plugin_system_name,
                          std::ostream      &output_stream,
                          const std::string &attachment_point)
      {
        // first output a graph node for the interface class as the central
        // hub of this plugin system, plotted as a square.
        //
        // we use the typeid name of the interface class to label
        // nodes within this plugin system, as they are unique among
        // all other plugin systems
        output_stream << std::string(typeid(InterfaceClass).name())
                      << " [label=\""
                      << plugin_system_name
                      << "\", height=.8,width=.8,shape=\"rect\",fillcolor=\"green\"]"
                      << std::endl;

        // then output the graph nodes for each plugin, with links to the
        // interface class and, as appropriate, from the SimulatorAccess class
        //
        // we would like to establish a predictable order of output here, but
        // plugins self-register via static global variables, and their
        // initialization order is not deterministic. consequently, let us
        // loop over all plugins first and put pointers to them into a
        // map with deterministic keys. as key, we use the declared name
        // of the plugin by which it is referred in the .prm file
        std::map<std::string, typename std::list<PluginInfo>::const_iterator>
        plugin_map;
        for (typename std::list<PluginInfo>::const_iterator p = plugins->begin();
             p != plugins->end(); ++p)
          plugin_map[std::get<0>(*p)] = p;

        // now output the information sorted by the plugin names
        for (typename std::map<std::string, typename std::list<PluginInfo>::const_iterator>::const_iterator
             p = plugin_map.begin();
             p != plugin_map.end(); ++p)
          {
            // take the name of the plugin and split it into strings of
            // 15 characters at most; then combine them
            // again using \n to make dot/neato show these parts of
            // the name on separate lines
            const std::vector<std::string> plugin_label_parts
              = Utilities::break_text_into_lines(p->first, 15);
            Assert (plugin_label_parts.size()>0, ExcInternalError());
            std::string plugin_name = plugin_label_parts[0];
            for (unsigned int i=1; i<plugin_label_parts.size(); ++i)
              plugin_name += "\\n" + plugin_label_parts[i];

            // next create a (symbolic) node name for this plugin. because
            // each plugin corresponds to a particular class, use the mangled
            // name of the class
            std::unique_ptr<InterfaceClass> instance (create_plugin (p->first, ""));
            const std::string node_name = typeid(*instance).name();

            // then output the whole shebang describing this node
            output_stream << node_name
                          << " [label=\""
                          << plugin_name
                          << "\", height=.8,width=.8,shape=\"circle\",fillcolor=\"lightblue\"];"
                          << std::endl;

            // next build connections from this plugin to the
            // interface class
            output_stream << node_name
                          << " -> "
                          << std::string(typeid(InterfaceClass).name())
                          << " [len=3, weight=50]"
                          << ';'
                          << std::endl;

            // finally see if this plugin is derived from
            // SimulatorAccess; if so, draw an arrow from SimulatorAccess
            // also to the plugin's name
            if (dynamic_cast<const SimulatorAccess<2>*>(instance.get()) != nullptr
                ||
                dynamic_cast<const SimulatorAccess<3>*>(instance.get()) != nullptr)
              output_stream << "SimulatorAccess"
                            << " -> "
                            << node_name
                            << " [style=\"dotted\", arrowhead=\"empty\", constraint=false, color=\"gray\", len=20, weight=0.1];"
                            << std::endl;
          }

        // as a last step, also draw a connection from the interface class
        // to the Simulator class, or whatever the calling function indicates
        // as the attachment point
        output_stream << std::string(typeid(InterfaceClass).name())
                      << " -> "
                      << attachment_point
                      << " [len=15, weight=50]"
                      << ';'
                      << std::endl;

        // end it with an empty line to make things easier to
        // read when looking over stuff visually
        output_stream << std::endl;
      }
    }
  }
}


#endif
