//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_table_h
#define __aspect__model_table_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A variable viscosity material model that reads the essential values of coefficients from
     * tables in input files.
     * 
     * This model is based on the paper Steinberger/Calderwood 2006:
     * "Models of large-scale viscous flow in the Earth's mantle with
     * contraints from mineral physics and surface observations"
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Steinberger: public MaterialModel::Interface<dim>
    {
      public:
        virtual double viscosity (const double temperature,
                                  const double pressure,
                                  const Point<dim> &position) const;

        virtual double reference_viscosity () const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const Point<dim> &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const;

        virtual bool is_compressible () const;
    };
  }
}

#endif
