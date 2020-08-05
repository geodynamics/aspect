/*
  Copyright (C) 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_boundary_velocity_residual_h
#define _aspect_postprocess_visualization_boundary_velocity_residual_h

#include <aspect/postprocess/visualization.h>
#include <aspect/boundary_velocity/gplates.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      using namespace dealii;
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the velocity residual between
       * the input velocity data and the velocity computed for the model domain.
       * The input velocity could either be an ascii data file or a file generated
       * by  the GPlates model.
       * This quantity only makes sense at the surface of the domain.
       * Thus, the value is set to zero in all the cells inside of the domain.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class BoundaryVelocityResidual
        : public DataPostprocessorVector<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          BoundaryVelocityResidual ();

          void initialize () override;

          /**
          * Evaluate the velocity residual for the current cell.
          *
          * @copydoc DataPostprocessorVector<dim>::evaluate_vector_field()
          */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const override;

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
           * Determines if the input ascii data file has velocity components in
           * spherical, i.e., (r, phi, theta) or in cartesian, i.e., (x, y, z)
           * coordinate system.
           */
          bool use_spherical_unit_vectors;

          /**
          * Determines if the model velocities are compared against ascii data
          * files, or a gplates model.
          */
          bool use_ascii_data;

          /**
          * Pointer to the gplates boundary velocity model
          */
          std::unique_ptr<BoundaryVelocity::internal::GPlatesLookup <dim> >  gplates_lookup;

          /**
           * Pointer to the ascii data
           */
          std::unique_ptr<Utilities::AsciiDataLookup<dim>> ascii_data_lookup;

          /**
           * Directory in which the input data files, i.e., GPlates model or ascii data
           * are present.
           */
          std::string data_directory;

          /**
           * Filename of the input Gplates model or ascii data file. For GPlates, the file names
           * can contain the specifiers %s and/or %c (in this order), meaning the name of the
           * boundary and the number of the data file time step.
           */
          std::string data_file_name;

          /**
           * Scale the input data by a scalar factor. Can be used to transform
           * the unit of the data (if they are not specified in SI units (m/s or
           * m/yr depending on the "Use years in output instead of seconds"
           * parameter).
           */
          double scale_factor;
      };
    }
  }
}

#endif