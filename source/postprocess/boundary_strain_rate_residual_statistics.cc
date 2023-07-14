/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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



#include <aspect/postprocess/boundary_strain_rate_residual_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    void
    BoundaryStrainRateResidualStatistics<dim>::initialize ()
    {
      strain_rate_data_lookup = std::make_unique<Utilities::StructuredDataLookup<dim>>(1, scale_factor);
      strain_rate_data_lookup->load_file(data_directory + data_file_name, this->get_mpi_communicator());
    }



    template <int dim>
    double
    BoundaryStrainRateResidualStatistics<dim>::get_data_surface_strain_rate (const Point<dim> &p) const
    {
      Point<dim> internal_position = p;

      if (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical)
        {
          const std::array<double,dim> spherical_position = this->get_geometry_model().
                                                            cartesian_to_natural_coordinates(p);

          for (unsigned int d = 0; d < dim; ++d)
            internal_position[d] = spherical_position[d];
        }

      const double data_surface_strain_rate = strain_rate_data_lookup->get_data(internal_position, 0);

      return data_surface_strain_rate;
    }



    template <int dim>
    std::pair<std::string,std::string>
    BoundaryStrainRateResidualStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula for the velocity.
      const Quadrature<dim-1> &quadrature_formula = this->introspection().face_quadratures.velocities;

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_gradients |
                                        update_JxW_values |
                                        update_quadrature_points);

      // This variable stores the strain rates computed at the quadrature points.
      std::vector<SymmetricTensor<2,dim>> strain_rate (fe_face_values.n_quadrature_points);

      std::map<types::boundary_id, double> local_max_eii;
      std::map<types::boundary_id, double> local_min_eii;
      std::map<types::boundary_id, double> local_eii_square_integral;
      std::map<types::boundary_id, double> local_boundary_area;

      std::set<types::boundary_id> boundary_indicators;
      boundary_indicators.insert(this->get_geometry_model().translate_symbolic_boundary_name_to_id("top"));

      AssertThrow (boundary_indicators.size() == 1,
                   ExcMessage("The <boundary strain rate residual statistics> plugin is only tested "
                              "to work correctly for a single boundary indicator, "
                              "but more than one was specified."));

      for (const auto boundary_id : boundary_indicators)
        {
          local_max_eii[boundary_id] = std::numeric_limits<double>::lowest();
          local_min_eii[boundary_id] = std::numeric_limits<double>::max();
        }

      // for every face that is part of the selected geometry boundary
      // and that is owned by this processor,
      // compute the maximum, minimum, and squared*area strain rate invariant residual
      // magnitude, and the face area.
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            for (const auto current_boundary_id : boundary_indicators)
              if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == current_boundary_id)
                {
                  fe_face_values.reinit (cell, f);

                  fe_face_values[this->introspection().extractors.velocities].get_function_symmetric_gradients
                  (this->get_solution(), strain_rate);

                  // determine the max, min, and area integral of squared strain rate residual on the face
                  // also determine the face area
                  double local_max = std::numeric_limits<double>::lowest();
                  double local_min = std::numeric_limits<double>::max();
                  double local_square_strain_rate = 0.0;
                  double local_fe_face_area = 0.0;
                  for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                    {
                      const Point<dim> point_at_surface = fe_face_values.quadrature_point(q);

                      // Extract the strain rate invariant output here.
                      const double e_ii = std::sqrt(std::fabs(second_invariant(deviator(strain_rate[q]))));

                      // Extract data strain rate invartiant.
                      double e_ii_data = get_data_surface_strain_rate(point_at_surface);

                      if (this->convert_output_to_years() == true)
                        e_ii_data = e_ii_data/year_in_seconds;

                      // The strain rate invariant residual is calculated here. Ignore the computation if we
                      // encounter nan or 1e300 in the input data.
                      double strain_invariant_residual = 0;
                      if (e_ii_data < 1e300 || !std::isnan(e_ii_data))
                        strain_invariant_residual = std::abs(e_ii_data - e_ii);
                      else
                        continue;

                      local_max = std::max(strain_invariant_residual,
                                           local_max);
                      local_min = std::min(strain_invariant_residual,
                                           local_min);
                      local_square_strain_rate += ((strain_invariant_residual * strain_invariant_residual) * fe_face_values.JxW(q));
                      local_fe_face_area += fe_face_values.JxW(q);
                    }
                  // then merge them with the min/max/squared strain rates
                  // and face areas we found for other faces with the same boundary indicator
                  local_max_eii[current_boundary_id] = std::max(local_max,
                                                                local_max_eii[current_boundary_id]);
                  local_min_eii[current_boundary_id] = std::min(local_min,
                                                                local_min_eii[current_boundary_id]);
                  local_eii_square_integral[current_boundary_id] += local_square_strain_rate;
                  local_boundary_area[current_boundary_id] += local_fe_face_area;
                }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_max_eii;
      std::map<types::boundary_id, double> global_min_eii;
      std::map<types::boundary_id, double> global_rms_eii;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        std::vector<double> local_max_values;
        std::vector<double> local_min_values;
        std::vector<double> local_eii_square_integral_values;
        std::vector<double> local_boundary_area_values;

        for (const auto p : boundary_indicators)
          {
            local_max_values.push_back (local_max_eii[p]);
            local_min_values.push_back (local_min_eii[p]);

            local_eii_square_integral_values.push_back (local_eii_square_integral[p]);
            local_boundary_area_values.push_back (local_boundary_area[p]);
          }

        // then collect contributions from all processors
        std::vector<double> global_max_values (local_max_values.size());
        Utilities::MPI::max (local_max_values, this->get_mpi_communicator(), global_max_values);
        std::vector<double> global_min_values (local_min_values.size());
        Utilities::MPI::min (local_min_values, this->get_mpi_communicator(), global_min_values);

        std::vector<double> global_eii_square_integral_values (local_eii_square_integral_values.size());
        Utilities::MPI::sum(local_eii_square_integral_values,this->get_mpi_communicator(),global_eii_square_integral_values);
        std::vector<double> global_boundary_area_values (local_boundary_area_values.size());
        Utilities::MPI::sum(local_boundary_area_values,this->get_mpi_communicator(),global_boundary_area_values);

        // and now take them apart into the global map again
        // At the same time, calculate the rms strain rate invariant for each boundary
        unsigned int index = 0;
        for (const auto boundary_id: boundary_indicators)
          {
            global_max_eii[boundary_id] = global_max_values[index];
            global_min_eii[boundary_id] = global_min_values[index];
            global_rms_eii[boundary_id] = std::sqrt(global_eii_square_integral_values[index] / global_boundary_area_values[index]);
            ++index;
          }
      }

      // now add the computed max, min, and rms strain rates to the statistics object
      // and create a single string that can be output to the screen
      const std::string units = (this->convert_output_to_years() == true) ? "1/yr" : "1/s";
      const double unit_scale_factor = (this->convert_output_to_years() == true) ? year_in_seconds : 1.0;
      std::ostringstream screen_text;
      unsigned int index = 0;

      for (std::map<types::boundary_id, double>::const_iterator
           max_eii = global_max_eii.begin(), min_eii = global_min_eii.begin(), rms = global_rms_eii.begin();
           max_eii != global_max_eii.end() && min_eii != global_min_eii.end() && rms != global_rms_eii.end();
           ++max_eii, ++min_eii, ++rms, ++index)
        {
          const std::string name_max = "Maximum strain rate invariant residual magnitude on boundary with indicator "
                                       + Utilities::int_to_string(max_eii->first)
                                       + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                     .translate_id_to_symbol_name (max_eii->first))
                                       + " (" + units + ")";
          statistics.add_value (name_max, max_eii->second*unit_scale_factor);
          const std::string name_min = "Minimum strain rate invariant residual magnitude on boundary with indicator "
                                       + Utilities::int_to_string(min_eii->first)
                                       + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                     .translate_id_to_symbol_name (min_eii->first))
                                       + " (" + units + ")";
          statistics.add_value (name_min, min_eii->second*unit_scale_factor);
          const std::string name_rms = "RMS strain rate invariant residual on boundary with indicator "
                                       + Utilities::int_to_string(rms->first)
                                       + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                     .translate_id_to_symbol_name (rms->first))
                                       + " (" + units + ")";
          statistics.add_value (name_rms, rms->second*unit_scale_factor);
          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name_max, 8);
          statistics.set_scientific (name_max, true);
          statistics.set_precision (name_min, 8);
          statistics.set_scientific (name_min, true);
          statistics.set_precision (name_rms, 8);
          statistics.set_scientific (name_rms, true);

          screen_text << max_eii->second *unit_scale_factor << " " << units << ", "
                      << min_eii->second *unit_scale_factor << " " << units << ", "
                      << rms->second *unit_scale_factor << " " << units
                      << (index == global_max_eii.size()-1 ? "" : ", ");

        }

      return std::pair<std::string, std::string> ("Max, min, and RMS residual strain rate invariant along boundary parts:",
                                                  screen_text.str());
    }



    template <int dim>
    void
    BoundaryStrainRateResidualStatistics<dim>::declare_parameters (ParameterHandler  &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Boundary strain rate residual statistics");
        {
          prm.declare_entry ("Data directory",
                             "$ASPECT_SOURCE_DIR/data/postprocess/boundary-strain-rate-residual/",
                             Patterns::DirectoryName (),
                             "The name of a directory that contains the ascii data. "
                             "This path may either be absolute (if starting with a "
                             "`/') or relative to the current directory. The path may also "
                             "include the special text `$ASPECT_SOURCE_DIR' which will be "
                             "interpreted as the path in which the ASPECT source files were "
                             "located when ASPECT was compiled. This interpretation allows, "
                             "for example, to reference files located in the `data/' subdirectory "
                             "of ASPECT.");
          prm.declare_entry ("Data file name", "box_3d_boundary_strain_rate.txt",
                             Patterns::Anything (),
                             "The file name of the input surface strain rate an ascii data. "
                             "The file has one column in addition to the coordinate columns "
                             "corresponding to the second invariant of strain rate. ");
          prm.declare_entry ("Scale factor", "1.",
                             Patterns::Double (),
                             "Scalar factor, which is applied to the model data. "
                             "You might want to use this to scale the input to a "
                             "reference model.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    BoundaryStrainRateResidualStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Boundary strain rate residual statistics");
        {
          // Get the path to the data files. If it contains a reference
          // to $ASPECT_SOURCE_DIR, replace it by what CMake has given us
          // as a #define
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          data_file_name    = prm.get ("Data file name");
          scale_factor      = prm.get_double ("Scale factor");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(BoundaryStrainRateResidualStatistics,
                                  "boundary strain rate residual statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the surface strain rate residual along the top boundary. The residual "
                                  "is the difference between the second invariant of the model strain rate and "
                                  "the second strain rate invariant read from the input data file. "
                                  "Currently, the strain residual statistics, i.e., min, max and the rms magnitude, "
                                  "are computed at the top surface.")
  }
}
