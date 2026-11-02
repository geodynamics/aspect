// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/postprocess/mobility_statistics.h>

#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    MobilityStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula for the velocity for the volume of cells
      const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.velocities;
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<Tensor<1,dim>> velocity_values(n_q_points);

      // create a quadrature formula for the velocity for the surface of cells
      const Quadrature<dim-1> &quadrature_formula_face
        = this->introspection().face_quadratures.velocities;
      const unsigned int n_q_points_face = quadrature_formula_face.size();

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula_face,
                                        update_values |
                                        update_JxW_values |
                                        update_quadrature_points);

      std::vector<Tensor<1,dim>> surface_velocity_values (n_q_points_face);

      double local_velocity_square_integral = 0;
      double local_velocity_square_integral_top_boundary = 0;
      double local_top_boundary_area = 0;

      const std::set<types::boundary_id>
      boundary_indicators
        = this->get_geometry_model().get_used_boundary_indicators ();
      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            // extract velocities
            fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                        velocity_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                local_velocity_square_integral += velocity_values[q].norm_square() *
                                                  fe_values.JxW(q);
              }

            // now for the surface cells only
            for (const unsigned int f : cell->face_indices())
              if (cell->face(f)->boundary_id() == top_boundary_id)
                {
                  fe_face_values.reinit (cell, f);
                  // extract velocities
                  fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                      surface_velocity_values);

                  // determine the squared velocity on the face
                  for (unsigned int q = 0; q < n_q_points_face; ++q)
                    {
                      const double JxW = fe_face_values.JxW(q);
                      local_velocity_square_integral_top_boundary += surface_velocity_values[q].norm_square() * JxW;
                      local_top_boundary_area += JxW;
                    }
                }
          }

      // now communicate to get the global values and convert to rms
      const double global_velocity_square_integral
        = Utilities::MPI::sum (local_velocity_square_integral, this->get_mpi_communicator());

      const double global_rms_vel = std::sqrt(global_velocity_square_integral) /
                                    std::sqrt(this->get_volume());

      const double global_top_velocity_square_integral = Utilities::MPI::sum(local_velocity_square_integral_top_boundary, this->get_mpi_communicator());
      const double global_top_boundary_area_integral = Utilities::MPI::sum(local_top_boundary_area, this->get_mpi_communicator());

      const double global_top_rms_vel = std::sqrt(global_top_velocity_square_integral) /
                                        std::sqrt(global_top_boundary_area_integral);

      // now add the computed rms velocities to the statistics object
      // and create a single string that can be output to the screen
      const double mobility  = global_top_rms_vel / global_rms_vel;
      const std::string name_mobility = "Mobility";

      statistics.add_value (name_mobility,
                            mobility);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name_mobility, 8);
      statistics.set_scientific (name_mobility, true);

      std::ostringstream output;
      output.precision(3);
      output << mobility;

      return std::pair<std::string, std::string> ("Mobility:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MobilityStatistics,
                                  "mobility statistics",
                                  "A postprocessor that computes some statistics about mobility "
                                  "following Tackley (2000) and Lourenco et al. (2020).")
  }
}
