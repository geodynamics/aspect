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
/*  $Id$  */


#ifndef __aspect__postprocess_visualization_h
#define __aspect__postprocess_visualization_h

#include <aspect/postprocess/interface.h>
#include <aspect/plugins.h>

#include <deal.II/base/thread_management.h>
#include <deal.II/numerics/data_postprocessor.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {

    }


    /**
     * A postprocessor that generates graphical output in periodic intervals
     * or every time step. The time interval between generating graphical
     * output is obtained from the parameter file.
     *
     * While this class acts as a plugin, i.e. as a postprocessor that can
     * be registered with the postprocessing manager class
     * aspect::Postprocess::Manager, this class at the same time also acts
     * as a manager for plugins itself, namely for classes derived from the
     * dealii::DataPostprocessor class that are used to output different
     * aspects of the solution, such as for example a computed seismic wave
     * speed from temperature, velocity and pressure.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Visualization : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Visualization ();

        /**
         * Destructor. Makes sure that any background thread that may still be
         * running writing data to disk finishes before the current object
         * is fully destroyed.
         */
        ~Visualization ();

        /**
         * Generate graphical output from the current solution.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Initialize this class for a given simulator. In addition to
         * calling the respective function from the base class, this
         * function also initializes all the visualization postprocessor
         * plugins.
         *
         * @param simulator A reference to the main simulator object to which the
         * postprocessor implemented in the derived class should be applied.
         **/
        virtual void initialize (const Simulator<dim> &simulator);


        /**
         * A function that is used to register visualization postprocessor objects
         * in such a way that the Manager can deal with all of them
         * without having to know them by name. This allows the files
         * in which individual postprocessors are implement to register
         * these postprocessors, rather than also having to modify the
         * Manage class by adding the new postprocessor class.
         *
         * @param name The name under which this visualization postprocessor is to
         * be called in parameter files.
         * @param description A text description of what this visualization plugin
         * does and that will be listed in the documentation of
         * the parameter file.
         * @param declare_parameters_function A pointer to a function
         * that declares the parameters for this postprocessor.
         * @param factory_function A pointer to a function that creates
         * such a postprocessor object and returns a pointer to it.
         **/
        static
        void
        register_visualization_postprocessor (const std::string &name,
                                              const std::string &description,
                                              void (*declare_parameters_function) (ParameterHandler &),
                                              DataPostprocessor<dim> *(*factory_function) ());

        /**
         * Declare the parameters this class takes through input files.
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

        /**
         * Save the state of this object.
         **/
        virtual
        void save (std::map<std::string, std::string> &status_strings) const;

        /**
         * Restore the state of the object.
         **/
        virtual
        void load (const std::map<std::string, std::string> &status_strings);

        /**
         * Serialize the contents of this class as far as they are not
         * read from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

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
         * Interval between the generation of graphical output. This
         * parameter is read from the input file and consequently is not part
         * of the state that needs to be saved and restored.
        *
        * For technical reasons, this value is stored as given in the
        * input file and upon use is either interpreted as seconds or
        * years, depending on how the global flag in the input parameter
        * file is set.
         */
        double output_interval;

        /**
         * A time (in years) after which the next time step should produce
         * graphical output again.
         */
        double next_output_time;

        /**
         * Consecutively counted number indicating the how-manyth time we
         * will create output the next time we get to it.
         */
        unsigned int output_file_number;

        /**
         * Graphical output format.
         */
        string output_format;

        /**
         * VTU file output supports grouping files from several CPUs
         * into one file using MPI I/O when writing on a parallel
         * filesystem. 0 means no grouping (and no parallel I/O).
         * 1 will generate one big file containing the whole solution.
         */
        unsigned int group_files;

        /**
         * Compute the next output time from the current one. In
         * the simplest case, this is simply the previous
         * next output time plus the interval, but in general
         * we'd like to ensure that it is larger than the current
         * time to avoid falling behind with next_output_time and
         * having to catch up once the time step becomes larger.
         */
        void set_next_output_time (const double current_time);

        /**
         * Handle to a thread that is used to write data in the
         * background. The background_writer() function runs
         * on this background thread.
         */
        Threads::Thread<void> background_thread;

        /**
         * A function that writes the text in the second argument
         * to a file with the name given in the first argument. The function
         * is run on a separate thread to allow computations to
         * continue even though writing data is still continuing.
         * The function takes over ownership of the arguments and deletes
        * them at the end of its work.
         */
        static
        void background_writer (const std::string *filename,
                                const std::string *file_contents);

        /**
         * A list of postprocessor objects that have been requested
         * in the parameter file.
         */
        std::list<std_cxx1x::shared_ptr<DataPostprocessor<dim> > > postprocessors;

    };
  }


  /**
   * Given a class name, a name, and a description for the parameter file for a postprocessor, register it with
   * the aspect::Postprocess::Manager class.
   *
   * @ingroup Postprocessing
   */
#define ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<dealii::DataPostprocessor<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::Postprocess::Manager<2>::register_visualization_postprocessor, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<dealii::DataPostprocessor<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Postprocess::Manager<3>::register_visualization_postprocessor, \
                                name, description); \
  }
}


#endif
