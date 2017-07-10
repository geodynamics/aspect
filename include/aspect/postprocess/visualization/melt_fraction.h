/*
  Copyright (C) 2013 - 2016 by the authors of the ASPECT code.

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

          virtual
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          /**
           * Parameters for anhydrous melting of peridotite after Katz, 2003
           */

          // for the solidus temperature
          double A1;   // °C
          double A2; // °C/Pa
          double A3; // °C/(Pa^2)

          // for the lherzolite liquidus temperature
          double B1;   // °C
          double B2;   // °C/Pa
          double B3; // °C/(Pa^2)

          // for the liquidus temperature
          double C1;   // °C
          double C2;  // °C/Pa
          double C3; // °C/(Pa^2)

          // for the reaction coefficient of pyroxene
          double r1;     // cpx/melt
          double r2;     // cpx/melt/GPa
          double M_cpx;  // mass fraction of pyroxenite

          // melt fraction exponent
          double beta;

          /**
           * Parameters for melting of pyroxenite after Sobolev et al., 2011
           */

          // for the melting temperature
          double D1;    // °C
          double D2;  // °C/Pa
          double D3; // °C/(Pa^2)

          // for the melt-fraction dependence of productivity
          double E1;
          double E2;
      };
    }
  }
}

#endif
