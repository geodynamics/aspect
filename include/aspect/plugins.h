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


#ifndef __aspect__plugins_h
#define __aspect__plugins_h


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <string>
#include <list>
#include <set>
#include <map>


namespace aspect
{
  namespace internal
  {
    /**
     * A namespace for the definition of classes that have to do with the
     * plugin architecture of Aspect.
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
        typedef
        std_cxx1x::tuple<std::string,
                  std::string,
                  void ( *) (ParameterHandler &),
                  InterfaceClass *( *) ()>
                  PluginInfo;

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
        if (plugins != 0)
          delete plugins;
        plugins = 0;
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
        if (plugins == 0)
          plugins = new std::list<PluginInfo>();

        // verify that the same name has not previously been
        // used to register a plugin, since we would then no
        // longer be able to identify the plugin
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          Assert (std_cxx1x::get<0>(*p) != name,
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
        Assert (plugins != 0,
                ExcMessage ("No plugins registered!?"));

        // get all names and put them into a data structure that keeps
        // them sorted
        std::set<std::string> names;
        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          names.insert (std_cxx1x::get<0>(*p));

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
          names_and_descriptions[std_cxx1x::get<0>(*p)] = std_cxx1x::get<1>(*p);;

        // then output it all
        typename std::map<std::string,std::string>::const_iterator
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
        Assert (plugins != 0,
                ExcMessage ("No postprocessors registered!?"));

        for (typename std::list<PluginInfo>::const_iterator
             p = plugins->begin();
             p != plugins->end(); ++p)
          (std_cxx1x::get<2>(*p))(prm);
      }



      template <typename InterfaceClass>
      InterfaceClass *
      PluginList<InterfaceClass>::
      create_plugin (const std::string  &name,
                     const std::string &documentation)
      {
        Assert (plugins != 0,
                ExcMessage ("No postprocessors registered!?"));
        Assert (name != "",
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
          if (std_cxx1x::get<0>(*p) == name)
            {
              InterfaceClass *i = std_cxx1x::get<3>(*p)();
              return i;
            }

        AssertThrow (false, ExcUnknownPlugin(name));
        return 0;
      }



      template <typename InterfaceClass>
      InterfaceClass *
      PluginList<InterfaceClass>::
      create_plugin (const std::string &name,
                     const std::string &documentation,
                     ParameterHandler &prm)
      {
        InterfaceClass *i = create_plugin(name, documentation);
        i->parse_parameters (prm);
        return i;
      }

    }
  }
}


#endif
