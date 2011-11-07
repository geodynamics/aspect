//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_simple_h
#define __aspect__model_simple_h

#include <aspect/model.h>

namespace aspect
{

  template <int dim>
  class MaterialModel_Simple: public MaterialModel<dim>
  {
    public:
      ~MaterialModel_Simple();

      virtual double eta (const double temperature, const double pressure, const dealii::Point<dim> &position) const;

      virtual double real_viscosity (const double                 temperature,
                                     const double                  pressure,
                                     const dealii::Point<dim> &position,
                                     const dealii::SymmetricTensor<2,dim> &strain_rate) const;

      virtual double specific_heat (const double temperature,
                                    const double pressure) const;

      virtual double density (const double temperature,
                              const double pressure,
                              const dealii::Point<dim> &position) const;

      virtual double compressibility (const double temperature,
                                      const double pressure,
                                      const dealii::Point<dim> &position) const;

      virtual double expansion_coefficient (const double temperature,
                                            const double pressure,
                                            const dealii::Point<dim> &position) const;
  };

}

#endif
