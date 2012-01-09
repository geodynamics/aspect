#ifndef __aspect__postprocess_TableVelocityStatistics_h
#define __aspect__postprocess_TableVelocityStatistics_h

#include <aspect/postprocess/interface.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the temperature.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class TableVelocityStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some temperature statistics for the Simple model.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
