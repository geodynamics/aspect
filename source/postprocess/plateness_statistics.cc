/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include <aspect/postprocess/plateness_statistics.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/mpi.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <utility>
#include <vector>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PlatenessStatistics<dim>::execute(TableHandler &statistics)
    {
      const Quadrature<dim-1> &quadrature_formula =
        this->introspection().face_quadratures.velocities;

      FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                       this->get_fe(),
                                       quadrature_formula,
                                       update_gradients |
                                       update_JxW_values |
                                       update_quadrature_points);

      std::vector<SymmetricTensor<2,dim>> strain_rate(fe_face_values.n_quadrature_points);

      const types::boundary_id top_boundary_id =
        this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      std::vector<std::pair<double,double>> local_deformation_area_pairs;
      double local_total_area = 0.0;
      double local_total_deformation = 0.0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->face(f)->at_boundary() &&
                cell->face(f)->boundary_id() == top_boundary_id)
              {
                fe_face_values.reinit(cell, f);

                fe_face_values[this->introspection().extractors.velocities]
                .get_function_symmetric_gradients(this->get_solution(), strain_rate);

                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const SymmetricTensor<2,dim> deviatoric_strain_rate =
                      Utilities::Tensors::consistent_deviator(strain_rate[q]);

                    const double e_ii =
                      std::sqrt(std::max(-Utilities::Tensors::consistent_second_invariant_of_deviatoric_tensor(deviatoric_strain_rate), 0.0));

                    const double area = fe_face_values.JxW(q);

                    local_deformation_area_pairs.emplace_back(e_ii, area);

                    local_total_area += area;
                    local_total_deformation += e_ii * area;
                  }
              }

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(this->get_mpi_communicator());

      const double global_total_area =
        Utilities::MPI::sum(local_total_area, this->get_mpi_communicator());
      const double global_total_deformation =
        Utilities::MPI::sum(local_total_deformation, this->get_mpi_communicator());

      AssertThrow(global_total_area > 0.0,
                  ExcMessage("The plateness statistics postprocessor found no nonzero area "
                             "on the top boundary. A nonzero surface area is required to "
                             "normalize F80 and F90."));
      AssertThrow(global_total_deformation > 0.0,
                  ExcMessage("The total deformation on the top boundary is zero. F80, F90, "
                             "and the corresponding plateness values are undefined when "
                             "there is no surface deformation."));

      const std::vector<std::vector<std::pair<double,double>>> gathered_deformation_area_pairs =
        Utilities::MPI::gather(this->get_mpi_communicator(), local_deformation_area_pairs);

      double f80 = 1.0;
      double f90 = 1.0;

      if (my_rank == 0)
        {
          std::vector<std::pair<double,double>> deformation_area_pairs;
          std::size_t n_deformation_area_pairs = 0;
          for (const auto &local_pairs : gathered_deformation_area_pairs)
            n_deformation_area_pairs += local_pairs.size();

          deformation_area_pairs.reserve(n_deformation_area_pairs);

          for (const auto &local_pairs : gathered_deformation_area_pairs)
            deformation_area_pairs.insert(deformation_area_pairs.end(),
                                          local_pairs.begin(),
                                          local_pairs.end());

          std::sort(deformation_area_pairs.begin(),
                    deformation_area_pairs.end(),
                    [](const auto &a, const auto &b)
          {
            return a.first > b.first;
          });

          double cumulative_deformation = 0.0;
          double cumulative_area = 0.0;
          bool found_f80 = false;
          bool found_f90 = false;

          for (const auto &entry : deformation_area_pairs)
            {
              cumulative_deformation += entry.first * entry.second;
              cumulative_area += entry.second;

              if (!found_f80 && cumulative_deformation >= 0.8 * global_total_deformation)
                {
                  f80 = cumulative_area / global_total_area;
                  found_f80 = true;
                }

              if (!found_f90 && cumulative_deformation >= 0.9 * global_total_deformation)
                {
                  f90 = cumulative_area / global_total_area;
                  found_f90 = true;
                  break;
                }
            }
        }

      MPI_Bcast(&f80, 1, MPI_DOUBLE, 0, this->get_mpi_communicator());
      MPI_Bcast(&f90, 1, MPI_DOUBLE, 0, this->get_mpi_communicator());

      const double p80 = 1.0 - f80 / reference_fraction;
      const double p90 = 1.0 - f90 / reference_fraction;

      statistics.add_value("F80 surface strain-rate invariant fraction", f80);
      statistics.add_value("F90 surface strain-rate invariant fraction", f90);
      statistics.add_value("Plateness p80", p80);
      statistics.add_value("Plateness p90", p90);

      statistics.set_precision("F80 surface strain-rate invariant fraction", 8);
      statistics.set_scientific("F80 surface strain-rate invariant fraction", false);
      statistics.set_precision("F90 surface strain-rate invariant fraction", 8);
      statistics.set_scientific("F90 surface strain-rate invariant fraction", false);
      statistics.set_precision("Plateness p80", 8);
      statistics.set_scientific("Plateness p80", false);
      statistics.set_precision("Plateness p90", 8);
      statistics.set_scientific("Plateness p90", false);

      std::ostringstream out;
      out << "F80=" << f80
          << ", F90=" << f90
          << ", p80=" << p80
          << ", p90=" << p90;

      return std::make_pair("Surface plateness statistics:", out.str());
    }


    template <int dim>
    void
    PlatenessStatistics<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Plateness statistics");
        {
          prm.declare_entry("Reference fraction", "0.6",
                            Patterns::Double(0.0, 1.0),
                            "Reference area fraction used to define plateness as "
                            "p = 1 - F/reference_fraction.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    PlatenessStatistics<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Plateness statistics");
        {
          reference_fraction = prm.get_double("Reference fraction");
          AssertThrow(reference_fraction > 0.0,
                      ExcMessage("The reference fraction must be greater than zero."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PlatenessStatistics,
                                  "plateness statistics",
                                  "A postprocessor that computes surface plateness diagnostics "
                                  "on the top boundary using the second invariant of the deviatoric "
                                  "strain-rate tensor. It computes F80 and F90, the fractional top-boundary "
                                  "area required to account for 80% and 90% of the total surface deformation, "
                                  "and also outputs p80 = 1 - F80/reference_fraction and "
                                  "p90 = 1 - F90/reference_fraction.")
  }
}
