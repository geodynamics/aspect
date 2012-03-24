#ifndef __aspect__initial_conditions_box_h
#define __aspect__initial_conditions_box_h

#include <aspect/initial_conditions/interface.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes a perturbed initial temperature field for
     * a box geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class PerturbedBox : public Interface<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature (const Point<dim> &position) const;
    };
  }
}

#endif
