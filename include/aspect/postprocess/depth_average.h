
//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__postprocess_visualization_h
#define __aspect__postprocess_visualization_h

#include <aspect/postprocess/interface.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that generates depth average output in periodic intervals
     * or every time step.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class DepthAverage : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        DepthAverage ();

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
         * Interval between the generation of output. This
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
         * Graphical output format.
         */
        string output_format;

        struct entry
        {
          double time, depth, value;
          template <class Archive>
          void serialize (Archive &ar, const unsigned int version);
        };

        /**
         * save all the past values
         */
        std::vector<entry> entries;

        /**
         * Compute the next output time from the current one. In
         * the simplest case, this is simply the previous
         * next output time plus the interval, but in general
         * we'd like to ensure that it is larger than the current
         * time to avoid falling behind with next_output_time and
         * having to catch up once the time step becomes larger.
         */
        void set_next_output_time (const double current_time);
    };
  }
}


#endif
