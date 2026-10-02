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


#ifndef _aspect_postprocess_visualization_melt_fraction_ball_h2o_h
#define _aspect_postprocess_visualization_melt_fraction_ball_h2o_h

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
      class MeltFractionBallH2O
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MeltFractionBallH2O ();

          /**
           * Melt fraction refers to the percentage of material that is molten for a given
           * @p temperature and @p pressure (assuming equilibrium conditions) for a given melting model.
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
           * Parameters for hydrous melting of peridotite after Katz, 2003.
           */

          // for the solidus temperature
          double A1; // °C
          double A2; // °C/Pa
          double A3; // °C/(Pa^2)

          // for the lherzolite liquidus temperature
          double B1; // °C
          double B2; // °C/Pa
          double B3; // °C/(Pa^2)

          // for the liquidus temperature
          double C1; // °C
          double C2; // °C/Pa
          double C3; // °C/(Pa^2)

          // for the reaction coefficient of pyroxene
          double r1; // cpx/melt
          double r2; // cpx/melt/GPa
          /**
           * This variable is read from the parameter file through a parameter called 'Mass fraction cpx'.
           */
          double M_cpx;

          // melt fraction exponent
          /**
           * The beta2 parameter is taken after Ball, 2022.
           * Equations 16, 17 and 18 from Katz et al 2003 together with their variables are defined here.
           */
          double beta1;
          double beta2;

          // eqn. 18 Katz et al 2003.
          double bulk_h2o_ppm;
          double D_H2O;
          double calc_X_H2O(const double F) const;

          // eqn. 16 Katz et al 2003.
          double K_H2O;
          double gamma_H2O;
          double calc_delta_T_H2O(const double F) const;

          // eqn. 17 Katz et al 2003. Water saturationn variables.
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

          // for the melting temperature
          double D1; // °C
          double D2; // °C/Pa
          double D3; // °C/(Pa^2)

          // for the melt-fraction dependence of productivity
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
