//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__gravity_model_vertical_h
#define __aspect__gravity_model_vertical_h

#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace GravityModel
  {
    using namespace dealii;

    /**
     * A class that describes gravity as a vector of constant
     * magnitude pointing vertically down.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class Vertical : public Interface<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const;
    };
  }
}

#endif
