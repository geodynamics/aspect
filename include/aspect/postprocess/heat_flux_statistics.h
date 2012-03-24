#ifndef __aspect__postprocess_heat_flux_statistics_h
#define __aspect__postprocess_heat_flux_statistics_h

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
    class HeatFluxStatistics : public Interface<dim>, public SimulatorAccess<dim>
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
