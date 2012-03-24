#ifndef __aspect__postprocess_table_velocity_statistics_h
#define __aspect__postprocess_table_velocity_statistics_h

#include <aspect/postprocess/interface.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the temperature.
     * This postprocessor can only be used with the MaterialModel::Table
     * model.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class TableVelocityStatistics : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some temperature statistics for the Table model.
         **/
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };
  }
}


#endif
