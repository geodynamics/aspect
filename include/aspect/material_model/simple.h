//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__model_simple_h
#define __aspect__model_simple_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that consists of globally constant values for all
     * material parameters except that the density decays linearly with the
     * temperature.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible. This is essentially
     * the material model used in the step-32 tutorial program.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Simple : public MaterialModel::Interface<dim>
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
        double eta;
	double thermal_expansion_coefficient;

					 /**
					  * The thermal conductivity.
					  */
	double k_value;
    };

  }
}

#endif
