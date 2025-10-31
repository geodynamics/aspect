/*
  Copyright (C) 2015 - 2018 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_seismic_melt_h
#define _aspect_postprocess_visualization_seismic_melt_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that outputs melt related
       * properties of the material model.
       */
      template <int dim>
      class MeltSeismicProperties
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          MeltSeismicProperties ();

          virtual
          std::vector<std::string>
          get_names () const;

          virtual
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const;

          virtual
          UpdateFlags
          get_needed_update_flags () const;

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
          std::vector<std::string> property_names;

          double solid_bulk_modulus_GPa;
          double liquid_bulk_modulus_GPa;
          double shear_modulus_GPa;
          double poisson;
          double p1;
          double p2;
          double p3;
          double p4;
          double p5;
          double p6;
      };

    }
  }
}

#endif
