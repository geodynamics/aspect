//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_simple_h
#define __aspect__model_simple_h

#include <aspect/model_base.h>

namespace aspect
{
  using namespace dealii;

  template <int dim>
  class MaterialModel_Simple: public MaterialModel<dim>
  {
    public:
      virtual double viscosity (const double temperature,
                                const double pressure,
                                const Point<dim> &position) const;

      virtual double specific_heat (const double temperature,
                                    const double pressure,
                                    const Point<dim> &position) const;

      virtual double density (const double temperature,
                              const double pressure,
                              const Point<dim> &position) const;

      virtual double compressibility (const double temperature,
                                      const double pressure,
                                      const Point<dim> &position) const;

      virtual bool is_compressible () const;

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
      double reference_density;
      double reference_temperature;
      double reference_eta;

  };

}

#endif
