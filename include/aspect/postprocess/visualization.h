//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__postprocess_visualization_h
#define __aspect__postprocess_visualization_h

#include <aspect/postprocess/interface.h>

#include <deal.II/base/thread_management.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that generates graphical output in periodic intervals
     * or every time step. The time interval between generating graphical
     * output is obtained from the parameter file.
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
         *
         * The communicator is a copy of the one that is used for
         * computations. However, it is not the same (but rather
         * a duplicate) since accessing it both from this thread
         * as well as from the main program has the potential to
         * produce race conditions.
         */
        static
        void background_writer (const std::string filename,
                                const std::string file_contents,
                                MPI_Comm          communicator);
    };
  }
}


#endif
