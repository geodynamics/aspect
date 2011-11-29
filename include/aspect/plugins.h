//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__plugins_h
#define __aspect__plugins_h


#include <deal.II/base/parameter_handler.h>
#include <string>


namespace aspect
{
  namespace internal
  {
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
    }
  }
}


#endif

class P;
