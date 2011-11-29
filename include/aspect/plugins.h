//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__plugins_h
#define __aspect__plugins_h


#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <string>
#include <list>


namespace aspect
{
  namespace internal
  {
				     /**
				      * A namespace for the definition of classes
				      * that have to do with the plugin architecture
				      * of Aspect.
				      */
    namespace Plugins
    {
      using namespace dealii;

      /**
       * An internal class that is used in the definition of the
       * ASPECT_REGISTER_* macros. Given a registration function, a classname,
       * a description of what it does, and a name for the parameter file,
       * it registers the model with the proper authorities.
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
         * Constructor. Given a pointer to a registration function
         * and name and description of the class, this
         * constructor registers the class passed as second
         * template argument.
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
         * A factory object that just creates object of the type registered
         * by this class.
         */
        static
        InterfaceClass *factory ()
        {
          return new ModelClass();
        }
      };


				       /**
					* A class that stores a list of registered
					* plugins for the given interface type.
					*/
      template <typename InterfaceClass>
      struct PluginList
      {
					   /**
					    * A type describing everything
					    * we need to know about a plugin.
					    *
					    * The entries in the tuple are:
					    * - The name by which it can be
					    *   selected.
					    * - A description of this plugin
					    *   that will show up in the
					    *   documentation in the
					    *   parameter file.
					    * - A function that can declare
					    *   the run-time parameters this
					    *   plugin takes from the parameter
					    *   file.
					    * - A function that can produce
					    *   objects of this plugin type.
					    */
	  typedef
	  std_cxx1x::tuple<std::string,
			   std::string,
			   void ( *) (ParameterHandler &),
			   InterfaceClass * ( *) ()>
	  PluginInfo;

					   /**
					    * A pointer to a list of
					    * all registered plugins.
					    *
					    * The object is a pointer
					    * rather than an object
					    * for the following
					    * reason: objects with
					    * static initializers
					    * (such as =0) are
					    * initialized before any
					    * objects for which one
					    * needs to run
					    * constructors.
					    * consequently, we can be
					    * sure that this pointer
					    * is set to zero before we
					    * ever try to register a
					    * postprocessor, and
					    * consequently whenever we
					    * run
					    * Manager::register_postprocessor,
					    * we need not worry
					    * whether we try to add
					    * something to this list
					    * before the lists's
					    * constructor has
					    * successfully run
					    */
	  static std::list<PluginInfo> *plugins;

					   /**
					    * Register a plugin by name, description,
					    * parameter declaration function, and
					    * factory function. See the discussion
					    * for the PluginInfo type above for
					    * more information on their meaning.
					    */
	  void register_plugin (const std::string &name,
				const std::string &description,
				void (*declare_parameters_function) (ParameterHandler &),
				InterfaceClass * (*factory_function) ())
	    {
					       // see if this is the first time we get into this
					       // function and if so initialize the variable above
	      if (plugins == 0)
		plugins = new std::list<PluginInfo>();

					       // now add one record to the list
	      plugins->push_back (PluginInfo(name,
					     description,
					     declare_parameters_function,
					     factory_function));
	    }
      };
    }
  }
}


#endif

class P;
