//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__postprocess_temperature_statistics_h
#define __aspect__postprocess_temperature_statistics_h

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
    class TemperatureStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        FILE *Tavout;
        /**
         * Evaluate the solution for some temperature statistics.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
