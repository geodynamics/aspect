/*
  Copyright (C) 2013 - 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_postprocess_visualization_melt_fraction_h
#define _aspect_postprocess_visualization_melt_fraction_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the melt fraction at every
       * point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class MeltFraction
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MeltFraction ();

          /**
           * Percentage of material that is molten for a given @p temperature and
           * @p pressure (assuming equilibrium conditions) for a given melting model.
           */
          double
          melt_fraction (const double temperature,
                         const double pressure,
                         const std::string &melting_model) const;

          /**
           * Evaluate the melt fraction for a given set of input data.
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Parameters for melting of peridotite after Katz, 2003.
           */

          /**
           * For the solidus temperature (respectively in °C, °C/Pa, °C/(Pa^2)).
           * These variables are read from the parameter file through parameters called 'A1', 'A2' and 'A3'.
           */
          double A1;
          double A2;
          double A3;

          /**
           * For the lherzolite liquidus temperature (respectively in °C, °C/Pa, °C/(Pa^2)).
           * These variables are read from the parameter file through parameters called 'B1', 'B2' and 'B3'.
           */
          double B1;
          double B2;
          double B3;

          /**
           * For the liquidus temperature (respectively in °C, °C/Pa, °C/(Pa^2)).
           * These variables are read from the parameter file through parameters called 'C1', 'C2' and 'C3'.
           */
          double C1;
          double C2;
          double C3;

          /**
           * For the reaction coefficient of pyroxene (units respectively cpx/melt and cpx/melt/GPa).
           * These variables are read from the parameter file through parameters called 'r1' and 'r2'.
           */
          double r1;
          double r2;

          /**
           * Mass fraction of pyroxenite.
           * This variable is read from the parameter file through a parameter called 'Mass fraction cpx'.
           */
          double M_cpx;

          /**
           * Equations 16, 17 and 18 from Katz et al 2003 together with their variables are defined here.
          */

          /**
           * Melt fraction exponents. The beta2 parameter is taken after Ball, 2022, the exponent when cpx is out.
           * These variables are read from the parameter file through parameters called 'beta1' and 'beta2'.
           */
          double beta1;
          double beta2;

          /**
           * These variables are read from the parameter file through parameters called 'bulk_h2o_ppm' and 'D_H2O'.
           * Eqn. 18 Katz et al (2003):
           *  Computes melt water concentration.
           *  Katz works with is %wt, so here ppms are converted to it.
           */
          double bulk_h2o_ppm;
          double D_H2O;
          double calc_X_H2O(const double F) const;

          /**
           * These variables are read from the parameter file through parameters called 'K_H2O' and 'gamma_H2O'.
           * Eqn. 16 Katz et al (2003):
           *  Computes the decrease of fusion temperature due to melt deluted water.
           */
          double K_H2O;
          double gamma_H2O;
          double calc_delta_T_H2O(const double F) const;

          /**
           * These variables are read from the parameter file through parameters called 'k1_H2O', 'k2_H2O' and 'lambda_H2O'.
           * Eqn. 17 Katz et al (2003):
           *  This function bounds the water concentration of melt to the water saturation, and recomputes melt fraction
           *  (using eqn. 19 Katz et al (2003)) with that value.
           */
          double k1_H2O;
          double k2_H2O;
          double lambda_H2O;
          double check_water_saturation(const double pressure,
                                        const double temperature,
                                        const double F,
                                        const double T_solidus,
                                        const double T_lherz_liquidus) const;

          /**
           * Parameters for melting of pyroxenite after Sobolev et al., 2011
           */

          /**
           * For the melting temperature (respectively in °C, °C/Pa, °C/(Pa^2)).
           * These variables are read from the parameter file through parameters called 'D1', 'D2' and 'D3'.
           */
          double D1;
          double D2;
          double D3;

          /**
           * For the melt-fraction dependence of productivity.
           * These variables are read from the parameter file through parameters called 'E1' and 'E2'.
           */
          double E1;
          double E2;

          /**
           * List of names of the melting models that are not peridotite.
           */
          std::vector<std::string> melting_model;
      };
    }
  }
}

#endif
