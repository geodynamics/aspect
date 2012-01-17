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
     * A material model that reads the essential values of coefficients from
     * tables in input files that describe their dependence as a function of
     * pressure and temperature.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Table: public MaterialModel::Interface<dim>
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

        virtual double thermal_diffusivity () const;

        virtual double density (const double temperature,
                                const double pressure,
                                const Point<dim> &position) const;

        virtual double Vp (const double temperature,
                           const double pressure) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const;

        virtual bool is_compressible () const;

        virtual double reference_density () const;

       double reference_gravity () const;

       double reference_thermal_alpha () const;

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
        double reference_rho;
        double reference_T;
        double reference_eta;
        double reference_kappa;
        double reference_alpha;
        double reference_g;
        std::string composition;
        std::string data_directory;

        /**
         * The thermal conductivity.
         */
      double k_value;
    };
  }
}

#endif
