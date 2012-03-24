#ifndef __aspect__gravity_model_radial_h
#define __aspect__gravity_model_radial_h

#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace GravityModel
  {
    using namespace dealii;

    /**
     * A class that describes gravity as a radial vector of constant
     * magnitude. The magnitude's value is read from the input file.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialConstant : public Interface<dim>
    {
      public:
        /**
         * Return the gravity vector as a function of position.
         */
        virtual Tensor<1,dim> gravity_vector (const Point<dim> &position) const;

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

      private:
        /**
         * Magnitude of the gravity vector.
         */
        double magnitude;
    };


    /**
     * A class that describes gravity as a radial vector with a magnitude that
     * is physically correct for an earth in which there are different densities
     * for the earth core and mantle. Specifically, at the core-mantle boundary,
     * gravity is assumed to be equal to 10.7 m/s^2 and it is 9.8 at the earth
     * surface; in between, it follows the behavior one would expect for a mantle
     * of constant density.
     *
     * This is the model used and discussed in the step-32 tutorial program.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class RadialEarthLike : public Interface<dim>
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
