/*
  Copyright (C) 2013 - 2021 by the authors of the ASPECT code.

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
           * Parameters for anhydrous melting of peridotite after Katz, 2003
           */

          // for the solidus temperature
          /**
           * This variable is read from the parameter file through a parameter called 'A1'.
           */
          double A1;   // °C
          /**
           * This variable is read from the parameter file through a parameter called 'A2'.
           */
          double A2; // °C/Pa
          /**
           * This variable is read from the parameter file through a parameter called 'A3'.
           */
          double A3; // °C/(Pa^2)

          // for the lherzolite liquidus temperature
          /**
           * This variable is read from the parameter file through a parameter called 'B1'.
           */
          double B1;   // °C
          /**
           * This variable is read from the parameter file through a parameter called 'B2'.
           */
          double B2;   // °C/Pa
          /**
           * This variable is read from the parameter file through a parameter called 'B3'.
           */
          double B3; // °C/(Pa^2)

          // for the liquidus temperature
          /**
           * This variable is read from the parameter file through a parameter called 'C1'.
           */
          double C1;   // °C
          /**
           * This variable is read from the parameter file through a parameter called 'C2'.
           */
          double C2;  // °C/Pa
          /**
           * This variable is read from the parameter file through a parameter called 'C3'.
           */
          double C3; // °C/(Pa^2)

          // for the reaction coefficient of pyroxene
          /**
           * This variable is read from the parameter file through a parameter called 'r1'.
           */
          double r1;     // cpx/melt
          /**
           * This variable is read from the parameter file through a parameter called 'r2'.
           */
          double r2;     // cpx/melt/GPa
          /**
           * This variable is read from the parameter file through a parameter called 'Mass fraction cpx'.
           */
          double M_cpx;  // mass fraction of pyroxenite

          // melt fraction exponent
          /**
           * This variable is read from the parameter file through a parameter called 'beta'.
           */
          double beta;

          /**
           * Parameters for melting of pyroxenite after Sobolev et al., 2011
           */

          // for the melting temperature
          /**
           * This variable is read from the parameter file through a parameter called 'D1'.
           */
          double D1;    // °C
          /**
           * This variable is read from the parameter file through a parameter called 'D2'.
           */
          double D2;  // °C/Pa
          /**
           * This variable is read from the parameter file through a parameter called 'D3'.
           */
          double D3; // °C/(Pa^2)

          // for the melt-fraction dependence of productivity
          /**
           * This variable is read from the parameter file through a parameter called 'E1'.
           */
          double E1;
          /**
           * This variable is read from the parameter file through a parameter called 'E2'.
           */
          double E2;
      };
    }
  }
}

#endif
