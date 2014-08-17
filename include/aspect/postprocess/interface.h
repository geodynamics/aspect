/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__postprocess_interface_h
#define __aspect__postprocess_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/std_cxx1x/shared_ptr.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/serialization/split_member.hpp>


namespace aspect
{
  using namespace dealii;

  template <int dim> class Simulator;
  template <int dim> class SimulatorAccess;


  /**
   * A namespace for everything to do with postprocessing solutions every time
   * step or every few time steps.
   *
   * @ingroup Postprocessing
   */
  namespace Postprocess
  {

    /**
     * This class declares the public interface of postprocessors.
     * Postprocessors must implement a function that can be called at the end
     * of each time step to evaluate the current solution, as well as
     * functions that save the state of the object and restore it (for
     * checkpoint/restart capabilities).
     *
     * Access to the data of the simulator is granted by the @p protected
     * member functions of the SimulatorAccess class, i.e., classes
     * implementing this interface will in general want to derive from both
     * this Interface class as well as from the SimulatorAccess class.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Does nothing but is virtual so that derived classes
         * destructors are also virtual.
         */
        virtual
        ~Interface ();

        /**
         * Initialize function.
         */
        virtual void initialize ();

        /**
         * Execute this postprocessor. Derived classes will implement this
         * function to do whatever they want to do to evaluate the solution at
         * the current time step.
         *
         * @param[in,out] statistics An object that contains statistics that
         * are collected throughout the simulation and that will be written to
         * an output file at the end of each time step. Postprocessors may
         * deposit data in these tables for later visualization or further
         * processing.
         *
         * @return A pair of strings that will be printed to the screen after
         * running the postprocessor in two columns; typically the first
         * column contains a description of what the data is and the second
         * contains a numerical value of this data. If there is nothing to
         * print, simply return two empty strings.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) = 0;

        /**
         * Declare the parameters this class takes through input files.
         * Derived classes should overload this function if they actually do
         * take parameters; this class declares a fall-back function that does
         * nothing, so that postprocessor classes that do not take any
         * parameters do not have to do anything at all.
         *
         * This function is static (and needs to be static in derived classes)
         * so that it can be called without creating actual objects (because
         * declaring parameters happens before we read the input file and thus
         * at a time when we don't even know yet which postprocessor objects
         * we need).
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation in this class does nothing, so that
         * derived classes that do not need any parameters do not need to
         * implement it.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


        /**
         * Save the state of this object to the argument given to this
         * function. This function is in support of checkpoint/restart
         * functionality.
         *
         * Derived classes can implement this function and should store their
         * state in a string that is deposited under a key in the map through
         * which the respective class can later find the status again when the
         * program is restarted. A legitimate key to store data under is
         * <code>typeid(*this).name()</code>. It is up to derived classes to
         * decide how they want to encode their state.
         *
         * The default implementation of this function does nothing, i.e., it
         * represents a stateless object for which nothing needs to be stored
         * at checkpoint time and nothing needs to be restored at restart
         * time.
         *
         * @param[in,out] status_strings The object into which implementations
         * in derived classes can place their status under a key that they can
         * use to retrieve the data.
         */
        virtual
        void save (std::map<std::string, std::string> &status_strings) const;

        /**
         * Restore the state of the object by looking up a description of the
         * state in the passed argument under the same key under which it was
         * previously stored.
         *
         * The default implementation does nothing.
         *
         * @param[in] status_strings The object from which the status will be
         * restored by looking up the value for a key specific to this derived
         * class.
         */
        virtual
        void load (const std::map<std::string, std::string> &status_strings);
    };






    /**
     * A class that manages all objects that provide functionality to
     * postprocess solutions. It declares run time parameters for input files,
     * reads their values from such an input file, manages a list of all
     * postprocessors selected in the input file, and upon request through the
     * execute() function calls them in turn.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Manager : public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Execute all of the postprocessor objects that have been requested
         * in the input file. These objects also fill the contents of the
         * statistics object.
         *
         * The function returns a concatenation of the text returned by the
         * individual postprocessors.
         */
        std::list<std::pair<std::string,std::string> >
        execute (TableHandler &statistics);

        /**
         * Go through the list of all postprocessors that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the desired type specified by the template
         * argument. If so, return a pointer to it. If no postprocessor is
         * active that matches the given type, return a NULL pointer.
         */
        template <typename PostprocessorType>
        PostprocessorType *
        find_postprocessor () const;

        /**
         * Declare the parameters of all known postprocessors, as well as of
         * ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which postprocessor objects will be created; then
         * let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Write the data of this object to a stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void save (Archive &ar,
                   const unsigned int version) const;

        /**
         * Read the data of this object from a stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void load (Archive &ar,
                   const unsigned int version);

        BOOST_SERIALIZATION_SPLIT_MEMBER()


        /**
         * A function that is used to register postprocessor objects in such a
         * way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * postprocessors are implement to register these postprocessors,
         * rather than also having to modify the Manager class by adding the
         * new postprocessor class.
         *
         * @param name The name under which this postprocessor is to be called
         * in parameter files.
         * @param description A text description of what this model does and
         * that will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that
         * declares the parameters for this postprocessor.
         * @param factory_function A pointer to a function that creates such a
         * postprocessor object and returns a pointer to it.
         */
        static
        void
        register_postprocessor (const std::string &name,
                                const std::string &description,
                                void (*declare_parameters_function) (ParameterHandler &),
                                Interface<dim> *(*factory_function) ());

        /**
         * Exception.
         */
        DeclException1 (ExcPostprocessorNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered postprocessors.");
      private:
        /**
         * A list of postprocessor objects that have been requested in the
         * parameter file.
         */
        std::list<std_cxx1x::shared_ptr<Interface<dim> > > postprocessors;
    };


    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <class Archive>
    void Manager<dim>::save (Archive &ar,
                             const unsigned int) const
    {
      // let all the postprocessors save their data in a map and then
      // serialize that
      std::map<std::string,std::string> saved_text;
      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        (*p)->save (saved_text);

      ar &saved_text;
    }


    template <int dim>
    template <class Archive>
    void Manager<dim>::load (Archive &ar,
                             const unsigned int)
    {
      // get the map back out of the stream; then let the postprocessors
      // that we currently have get their data from there. note that this
      // may not be the same set of postprocessors we had when we saved
      // their data
      std::map<std::string,std::string> saved_text;
      ar &saved_text;

      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        (*p)->load (saved_text);
    }

    /**
     * Go through the list of all postprocessors that have been selected in
     * the input file (and are consequently currently active) and see if one
     * of them has the desired type specified by the template argument. If so,
     * return a pointer to it. If no postprocessor is active that matches the
     * given type, return a NULL pointer.
     */
    template <int dim>
    template <typename PostprocessorType>
    inline
    PostprocessorType *
    Manager<dim>::find_postprocessor () const
    {
      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        if (PostprocessorType *x = dynamic_cast<PostprocessorType *> ( (*p).get()) )
          return x;
      return NULL;
    }


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a postprocessor, register it with the aspect::Postprocess::Manager
     * class.
     *
     * @ingroup Postprocessing
     */
#define ASPECT_REGISTER_POSTPROCESSOR(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_POSTPROCESSOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Postprocess::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::Postprocess::Manager<2>::register_postprocessor, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Postprocess::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Postprocess::Manager<3>::register_postprocessor, \
                                name, description); \
  }
  }
}


#endif
