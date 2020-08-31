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

#include <aspect/simulator.h>
#include <aspect/postprocess/visualization/boundary_velocity_residual.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      BoundaryVelocityResidual<dim>::
      BoundaryVelocityResidual ()
        :
        DataPostprocessorVector<dim> ("boundary_velocity_residual",
                                      update_values | update_quadrature_points | update_gradients)
      {}



      template <int dim>
      void
      BoundaryVelocityResidual<dim>::initialize ()
      {

        if (use_ascii_data)
          {
            // The input ascii table contains dim data columns (velocity components) in addition to the coordinate columns.
            ascii_data_lookup = std_cxx14::make_unique<Utilities::AsciiDataLookup<dim> >(dim, scale_factor);
            ascii_data_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
          }
        else
          {
            // The variables pointone and pointtwo are used in gplates to find the 2D plane in which the model lies.
            Tensor<1,2> pointone;
            Tensor<1,2> pointtwo;

            // These values are not used for 3D geometries and are thus set to the default.
            pointone[0] = numbers::PI/2;
            pointone[1] = 0;
            pointtwo[0] = numbers::PI/2;
            pointtwo[1] = numbers::PI/2;

            gplates_lookup =  std_cxx14::make_unique<BoundaryVelocity::internal::GPlatesLookup<dim> > (
                                pointone, pointtwo);

            gplates_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
          }
      }



      template <int dim>
      void
      BoundaryVelocityResidual<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        Assert ((computed_quantities[0].size() == dim), ExcInternalError());
        auto cell = input_data.template get_cell<DoFHandler<dim> >();

        const double velocity_scaling_factor =
          this->convert_output_to_years() ? year_in_seconds : 1.0;

        for (unsigned int q=0; q<computed_quantities.size(); ++q)
          for (unsigned int d = 0; d < dim; ++d)
            computed_quantities[q](d)= 0.;

        // We only want the output at the top boundary, so only compute it if the current cell
        // has a face at the top boundary.
        bool cell_at_top_boundary = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f) &&
              (this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top"))
            cell_at_top_boundary = true;


        if (cell_at_top_boundary)
          {
            for (unsigned int q=0; q<computed_quantities.size(); ++q)
              {
                Tensor<1,dim> data_velocity;
                if (use_ascii_data)
                  {
                    Point<dim> internal_position = input_data.evaluation_points[q];

                    if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
                      {
                        const std::array<double,dim> spherical_position = this->get_geometry_model().
                                                                          cartesian_to_natural_coordinates(input_data.evaluation_points[q]);

                        for (unsigned int d = 0; d < dim; ++d)
                          internal_position[d] = spherical_position[d];
                      }

                    for (unsigned int d = 0; d < dim; ++d)
                      data_velocity[d] = ascii_data_lookup->get_data(internal_position, d);

                    if (use_spherical_unit_vectors == true)
                      data_velocity = Utilities::Coordinates::spherical_to_cartesian_vector(data_velocity, input_data.evaluation_points[q]);
                  }

                else
                  {
                    data_velocity =  gplates_lookup->surface_velocity(input_data.evaluation_points[q]);
                  }

                for (unsigned int d = 0; d < dim; ++d)
                  computed_quantities[q](d) = data_velocity[d] - input_data.solution_values[q][d] * velocity_scaling_factor;
              }
          }
      }



      template <int dim>
      void
      BoundaryVelocityResidual<dim>::declare_parameters (ParameterHandler  &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Boundary velocity residual");
            {
              prm.declare_entry ("Data directory",
                                 "$ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/",
                                 Patterns::DirectoryName (),
                                 "The name of a directory that contains the GPlates model or the "
                                 "ascii data. This path may either be absolute (if starting with a "
                                 "`/') or relative to the current directory. The path may also "
                                 "include the special text `$ASPECT_SOURCE_DIR' which will be "
                                 "interpreted as the path in which the ASPECT source files were "
                                 "located when ASPECT was compiled. This interpretation allows, "
                                 "for example, to reference files located in the `data/' subdirectory "
                                 "of ASPECT.");
              prm.declare_entry ("Data file name", "current_day.gpml",
                                 Patterns::Anything (),
                                 "The file name of the input velocity as a GPlates model or an ascii data. "
                                 "For the GPlates model, provide file in the same format as described "
                                 "in the 'gplates' boundary velocity plugin. "
                                 "For the ascii data, provide file in the same format as described in "
                                 " 'ascii data' initial composition plugin." );
              prm.declare_entry ("Scale factor", "1.",
                                 Patterns::Double (),
                                 "Scalar factor, which is applied to the model data. "
                                 "You might want to use this to scale the input to a "
                                 "reference model. Another way to use this factor is to "
                                 "convert units of the input files. For instance, if you "
                                 "provide velocities in cm/yr set this factor to 0.01.");
              prm.declare_entry ("Use spherical unit vectors", "false",
                                 Patterns::Bool (),
                                 "Specify velocity as r, phi, and theta components "
                                 "instead of x, y, and z. Positive velocities point up, east, "
                                 "and north (in 3D) or out and clockwise (in 2D). "
                                 "This setting only makes sense for spherical geometries."
                                 "GPlates data is always interpreted to be in east and north directions "
                                 "and is not affected by this parameter.");
              prm.declare_entry ("Use ascii data", "false",
                                 Patterns::Bool (),
                                 "Use ascii data files (e.g., GPS) for computing residual velocities "
                                 "instead of GPlates data.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      BoundaryVelocityResidual<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        prm.enter_subsection("Visualization");
        {
          prm.enter_subsection("Boundary velocity residual");
          {
            // Get the path to the data files. If it contains a reference
            // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
            // as a #define
            data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
            data_file_name    = prm.get ("Data file name");
            scale_factor      = prm.get_double ("Scale factor");

            use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
            use_ascii_data = prm.get_bool("Use ascii data");

            if (use_spherical_unit_vectors)
              AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                           ExcMessage ("Spherical unit vectors should not be used "
                                       "when geometry model is not spherical."));
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
        prm.leave_subsection();
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(BoundaryVelocityResidual,
                                                  "boundary velocity residual",
                                                  "A visualization output object that generates output for the velocity "
                                                  "residual at the top surface. The residual is computed at each point at the "
                                                  "surface as the difference between the modeled velocities and the input "
                                                  "data velocities for each vector component. The user has an option to choose "
                                                  "the input data as ascii data files (e.g. GPS velocities) with columns "
                                                  "in the same format as described for the 'ascii data' initial temperature plugin "
                                                  "or a velocity field computed from the GPlates program as described in the gplates "
                                                  "boundary velocity plugin. ")
    }
  }
}