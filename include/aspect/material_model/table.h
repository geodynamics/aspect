//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
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
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double temperature,
                                  const double pressure,
                                  const Point<dim> &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        /**
         * A reference thermal diffusivity $\kappa$. $\kappa$ is related to the thermal
         * conductivity $k$ as $\kappa = k/(rho c_p)$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        double reference_thermal_diffusivity () const;

        /**
         * A reference thermal expansion coefficient $\alpha$.
         *
         * The value here is not used in the computation of things but only in
         * postprocessing the solution when we want dimension-less
         * quantities.
         */
        double reference_thermal_alpha () const;

        /**
        * A reference thermal expansion coefficient $\alpha$.
        *
        * The value here is not used in the computation of things but only in
        * postprocessing the solution when we want dimension-less
        * quantities.
        */
        double reference_cp () const;
        /**
         * @}
         */

        /**
         * @name Auxiliary material properties used for postprocessing
         * @{
         */
        virtual double seismic_Vp (const double temperature,
                                   const double pressure) const;

        virtual double seismic_Vs (const double temperature,
                                   const double pressure) const;

        virtual std::string datadir () const;

        virtual unsigned int thermodynamic_phase (const double temperature,
                                                  const double pressure) const;


        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
        /**
         * @}
         */

      private:
        double reference_rho;
        double reference_T;
        double reference_kappa;
        double reference_specific_heat;
        double reference_alpha;
        std::string composition;
        std::string data_directory;
        bool ComputePhases;
        bool Compressible;

        std::string ViscosityModel;
        double reference_eta;
        double ExponentialT;
        double ExponentialP;
        /**
         * The thermal conductivity.
         */
        double k_value;
    };
  }
}

#endif
