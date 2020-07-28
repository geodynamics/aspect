/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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
      // The input gpaltes file with default values of pointone and pointtwo
	    pointone[0] = 1.570796; pointone[1] = 0;
	    pointtwo[0] = 1.570796; pointtwo[1] = 1.570796;

	   gplates_lookup =  std_cxx14::make_unique<BoundaryVelocity::internal::GPlatesLookup<dim> > (
			   pointone, pointtwo);

	   gplates_lookup->load_file(gplates_data_directory + gplates_data_file_name, this->get_mpi_communicator());

	 // The input ascii table contains one data column (velocity components) in addition to the coordinate columns.
	   ascii_data_lookup = std_cxx14::make_unique<Utilities::AsciiDataLookup<dim> >(dim, ascii_scale_factor);
	   ascii_data_lookup->load_file(ascii_data_directory + ascii_data_file_name, this->get_mpi_communicator());
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
        	computed_quantities[q](d)= 0;

        // We only want to output dynamic topography at the top and bottom
        // boundary, so only compute it if the current cell has
        // a face at the top or bottom boundary.
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
		    if (use_spherical_unit_vectors && use_gps_data)
			{
		      const std::array<double,dim> cartesian_position = this->get_geometry_model().
		    		  	  	  	  cartesian_to_natural_coordinates(input_data.evaluation_points[q]);

		      Point<dim> internal_position;

		      for (unsigned int d = 0; d < dim; ++d)
		    	internal_position[d] = cartesian_position[d];

		      for (unsigned int d = 0; d < dim; ++d)
		    	data_velocity[d] = ascii_data_lookup->get_data(internal_position, d);

		      data_velocity = Utilities::Coordinates::spherical_to_cartesian_vector(data_velocity, input_data.evaluation_points[q]);

		      for (unsigned int d = 0; d < dim; ++d)
		        computed_quantities[q](d) = (data_velocity[d] - input_data.solution_values[q][d] * velocity_scaling_factor) ;
			}
		    else if (use_gps_data && !use_spherical_unit_vectors)
		    {
		      for (unsigned int d = 0; d < dim; ++d)
		      {
			    data_velocity[d] = ascii_data_lookup->get_data(input_data.evaluation_points[q], d);
			    computed_quantities[q](d) = (data_velocity[d] - input_data.solution_values[q][d] * velocity_scaling_factor) ;
		      }
		    }
		    else
		    {
		      data_velocity =  gplates_lookup->surface_velocity(input_data.evaluation_points[q]);
		      for (unsigned int d = 0; d < dim; ++d)
		        computed_quantities[q](d) = (data_velocity[d] - input_data.solution_values[q][d]) * velocity_scaling_factor;
		    }
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
		    prm.enter_subsection("Surface properties");
		    {
			prm.declare_entry ("Ascii data directory",
					           "$ASPECT_SOURCE_DIR//data/boundary-velocity/ascii-data/test/",
							   Patterns::DirectoryName (),
							   "The name of a directory that contains the model data. This path "
							   "may either be absolute (if starting with a `/') or relative to "
							   "the current directory. The path may also include the special "
							   "text `$ASPECT_SOURCE_DIR' which will be interpreted as the path "
							   "in which the ASPECT source files were located when ASPECT was "
							   "compiled. This interpretation allows, for example, to reference "
							   "files located in the `data/' subdirectory of ASPECT. ");
			prm.declare_entry ("Ascii data file name",
							   "shell_3d_top_spherical.0.txt",
							   Patterns::Anything (),
							   "The file name of the model data. Provide file in format: "
							   "(Velocity file name).\\%s\\%d where \\%s is a string specifying "
							   "the boundary of the model according to the names of the boundary "
							   "indicators (of the chosen geometry model).\\%d is any sprintf integer "
							   "qualifier, specifying the format of the current file number. ");
			prm.declare_entry ("Scale factor", "1.",
							   Patterns::Double (),
							   "Scalar factor, which is applied to the model data. "
							   "You might want to use this to scale the input to a "
							   "reference model. Another way to use this factor is to "
							   "convert units of the input files. For instance, if you "
							   "provide velocities in cm/yr set this factor to 0.01.");
			prm.declare_entry ("GPlates data directory",
							   "$ASPECT_SOURCE_DIR/data/boundary-velocity/gplates/",
							   Patterns::DirectoryName (),
							   "The name of a directory that contains the model data. This path "
							   "may either be absolute (if starting with a `/') or relative to "
							   "the current directory. The path may also include the special "
							   "text `$ASPECT_SOURCE_DIR' which will be interpreted as the path "
							   "in which the ASPECT source files were located when ASPECT was "
							   "compiled. This interpretation allows, for example, to reference "
							   "files located in the `data/' subdirectory of ASPECT. ");
			prm.declare_entry ("GPlates data file name", "current_day.gpml",
							   Patterns::Anything (),
							   "The file name of the model data. Provide file in format: "
							   "(Velocity file name).\\%s\\%d where \\%s is a string specifying "
							   "the boundary of the model according to the names of the boundary "
							   "indicators (of the chosen geometry model).\\%d is any sprintf integer "
							   "qualifier, specifying the format of the current file number. ");
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
							   "This setting only makes sense for spherical geometries.");
		    prm.declare_entry ("Use GPS data", "false",
							   Patterns::Bool (),
							   "Use GPS data for residual instead of the GPlates velocity data.");
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
			prm.enter_subsection("Surface properties");
			{
			// Get the path to the data files. If it contains a reference
			// to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
			// as a #define
			ascii_data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Ascii data directory"));
			ascii_data_file_name    = prm.get ("Ascii data file name");
			ascii_scale_factor      = prm.get_double ("Scale factor");

			gplates_data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("GPlates data directory"));
			gplates_data_file_name    = prm.get ("GPlates data file name");
			gplates_scale_factor      = prm.get_double ("Scale factor");

			use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
			use_gps_data = prm.get_bool("Use GPS data");

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
												  "A visualization output object that generates output for the velocity"
												  "residual at the top surface. The residual is computed at each point as "
												  "the difference of the modeled velocities from the input velocities at"
												  "the surface for each vector component. The user has an option to choose"
												  "the input data as the GPS velocities in the ASCII data format with columns "
												  "in the same format as described for the initial temperature or compositions"
												  "or a velocity filed computed from the GPlates program. ")
    }
  }
}
