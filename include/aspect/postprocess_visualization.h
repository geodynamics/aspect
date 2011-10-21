//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__postprocess_visualization_h
#define __aspect__postprocess_visualization_h

#include <aspect/postprocess_base.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that generates graphical output in periodic intervals
     * or every time step.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Visualization : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        Visualization ();

        /**
         * Generate graphical output from the current solution.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);


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

      private:
        /**
         * Interval in years between the generation of graphical output. This
         * parameter is read from the input file and consequently is not part
         * of the state that needs to be saved and restored.
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
    };
  }
}


#endif
