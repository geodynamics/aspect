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
  using namespace dealii;
  
  /**
   * A base class for parameterizations of material models.Classes derived from
   * this class will need to implement functions that provide material parameters
   * such as the viscosity, density, etc, typically as a function of position,
   * temperature and pressure at that location.
   */
  template <int dim>
  class MaterialModel
  {
    public:
      /**
       * Destructor. Made virtual to enforce that derived classes also have
       * virtual destructors.
       */
      virtual ~MaterialModel();

      /**
       * Return the viscosity $\eta$ of the model as a function of temperature,
       * pressure and position.
       */
      virtual double viscosity (const double      temperature, 
				const double      pressure, 
				const Point<dim> &position) const = 0;

      virtual double real_viscosity (const double      temperature,
                                     const double      pressure,
                                     const Point<dim> &position,
                                     const SymmetricTensor<2,dim> &strain_rate) const = 0;

      /**
       * Return the specific heat (i.e. $c_P$) of the model as a function of temperature,
       * pressure and position.
       */
      virtual double specific_heat (const double      temperature,
                                    const double      pressure,
				    const Point<dim> &position) const = 0;

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

      static
      MaterialModel<dim> *
      create( const std::string name);

  };

}

#endif
