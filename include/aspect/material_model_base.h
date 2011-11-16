//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_h
#define __aspect__model_h

#include <deal.II/numerics/vectors.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with modeling
   * convecting material, including descriptions of material parameters such
   * as viscosities, densities, etc.
   */
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A base class for parameterizations of material models. Classes derived from
     * this class will need to implement functions that provide material parameters
     * such as the viscosity, density, etc, typically as a function of position,
     * temperature and pressure at that location.
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Return the viscosity $\eta$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double viscosity (const double      temperature,
                                  const double      pressure,
                                  const Point<dim> &position) const = 0;

        /**
         * Return a reference value typical of the viscosities that
         * appear in this model. This value is not actually used in the
         * material description itself, but is used in scaling variables
         * to the same numerical order of magnitude when solving linear
         * systems. Specifically, the reference viscosity appears in
         * the factor scaling the pressure against the velocity.
         */
        virtual double reference_viscosity () const = 0;

        /**
         * Return the specific heat $c_P$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double specific_heat (const double      temperature,
                                      const double      pressure,
                                      const Point<dim> &position) const = 0;

        /**
         * Return the thermal conductivity $k$ of the model.
         */
        virtual double thermal_conductivity () const = 0;

        /**
         * Return the density $\rho$ of the model as a function of temperature,
         * pressure and position.
         */
        virtual double density (const double      temperature,
                                const double      pressure,
                                const Point<dim> &position) const = 0;

        /**
         * Return the compressibility coefficient of the model as a function of temperature,
         * pressure and position.
         */
        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const = 0;

        /**
         * Return whether the model is compressible or not.
         **/
        virtual bool is_compressible () const = 0;

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
    };


    /**
     * A function that given the name of a model returns a pointer to an object
     * that describes it. Ownership of the pointer is transferred to the caller.
     */
    template <int dim>
    Interface<dim> *
    create (const std::string &name);
  }

}

#endif
