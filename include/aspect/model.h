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

  template <int dim>
  class MaterialModel
  {
    public:
      virtual ~MaterialModel();

      virtual double eta (const double temperature, const double pressure, const dealii::Point<dim> &position) const = 0;

      virtual double real_viscosity (const double                 temperature,
                                     const double                  pressure,
                                     const dealii::Point<dim> &position,
                                     const dealii::SymmetricTensor<2,dim> &strain_rate) const = 0;

      // aka rho-cp:
      virtual double specific_heat (const double temperature,
                                    const double pressure) const = 0;

      virtual double density (const double temperature,
                              const double pressure,
                              const dealii::Point<dim> &position) const = 0;

      virtual double compressibility (const double temperature,
                                      const double pressure,
                                      const dealii::Point<dim> &position) const = 0;

      virtual double expansion_coefficient (const double temperature,
                                            const double pressure,
                                            const dealii::Point<dim> &position) const = 0;

      /**
       * Declare the parameters this class takes through input files.
       */
      static
      void
      declare_parameters (dealii::ParameterHandler &prm);

      /**
       * Read the parameters this class declares from the parameter
       * file.
       */
      virtual
      void
      parse_parameters (dealii::ParameterHandler &prm);

  };

}

#endif
