//-------------------------------------------------------------
//    $Id: heat_flux_statistics.h 573 2012-01-03 08:25:46Z geenen $
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__postprocess_table_heat_flux_statistics_h
#define __aspect__postprocess_table_heat_flux_statistics_h

#include <aspect/postprocess/interface.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the heat_flux.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class TableHeatfluxStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some heat_flux statistics.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
