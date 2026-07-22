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
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      /**
       * Compute the surface-area fractions that contain 80% and 90% of the
       * total deformation without collecting the quadrature-point data on one
       * MPI process. The algorithm bisects the global range of strain-rate
       * invariants and uses global sums to determine how much deformation lies
       * above each trial threshold. At the final threshold, it includes only
       * the fraction of the equal-strain-rate area needed to reach the target
       * deformation exactly.
       *
       * Each entry in @p local_deformation_area_pairs contains a strain-rate
       * invariant as its first value and its quadrature-point area as its
       * second value.
       */
      std::array<double,2>
      compute_deformation_area_fractions(const std::vector<std::pair<double,double>> &local_deformation_area_pairs,
                                         const double global_total_deformation,
                                         const double global_total_area,
                                         const MPI_Comm mpi_communicator)
      {
        const std::array<double,2> target_fractions = {{0.8, 0.9}};

        double local_minimum_strain_rate = std::numeric_limits<double>::max();
        double local_maximum_strain_rate = std::numeric_limits<double>::lowest();
        for (const auto &entry : local_deformation_area_pairs)
          {
            local_minimum_strain_rate = std::min(local_minimum_strain_rate, entry.first);
            local_maximum_strain_rate = std::max(local_maximum_strain_rate, entry.first);
          }

        const double global_minimum_strain_rate =
          Utilities::MPI::min(local_minimum_strain_rate, mpi_communicator);
        const double global_maximum_strain_rate =
          Utilities::MPI::max(local_maximum_strain_rate, mpi_communicator);

        if (global_minimum_strain_rate == global_maximum_strain_rate)
          return target_fractions;

        std::array<double,2> lower_thresholds;
        std::array<double,2> upper_thresholds;
        lower_thresholds.fill(std::nextafter(global_minimum_strain_rate,
                                             std::numeric_limits<double>::lowest()));
        upper_thresholds.fill(std::nextafter(global_maximum_strain_rate,
                                             std::numeric_limits<double>::max()));

        for (unsigned int iteration=0; iteration<std::numeric_limits<double>::digits; ++iteration)
          {
            double test_thresholds[2];
            for (unsigned int i=0; i<target_fractions.size(); ++i)
              {
                test_thresholds[i] =
                  (lower_thresholds[i] > 0.0
                   ? std::sqrt(lower_thresholds[i] * upper_thresholds[i])
                   : 0.5 * (lower_thresholds[i] + upper_thresholds[i]));

                if (test_thresholds[i] <= lower_thresholds[i] ||
                    test_thresholds[i] >= upper_thresholds[i])
                  test_thresholds[i] = lower_thresholds[i];
              }

            double local_deformation_above_threshold[2] = {0.0, 0.0};
            for (const auto &entry : local_deformation_area_pairs)
              for (unsigned int i=0; i<target_fractions.size(); ++i)
                if (entry.first > test_thresholds[i])
                  local_deformation_above_threshold[i] += entry.first * entry.second;

            double global_deformation_above_threshold[2];
            Utilities::MPI::sum(local_deformation_above_threshold,
                                mpi_communicator,
                                global_deformation_above_threshold);

            for (unsigned int i=0; i<target_fractions.size(); ++i)
              if (global_deformation_above_threshold[i] > target_fractions[i] * global_total_deformation)
                lower_thresholds[i] = test_thresholds[i];
              else
                upper_thresholds[i] = test_thresholds[i];
          }

        double local_cutoff_strain_rates[2] = {std::numeric_limits<double>::max(),
                                               std::numeric_limits<double>::max()
                                              };
        for (const auto &entry : local_deformation_area_pairs)
          for (unsigned int i=0; i<target_fractions.size(); ++i)
            if (entry.first > lower_thresholds[i])
              local_cutoff_strain_rates[i] = std::min(local_cutoff_strain_rates[i], entry.first);

        double cutoff_strain_rates[2];
        Utilities::MPI::min(local_cutoff_strain_rates,
                            mpi_communicator,
                            cutoff_strain_rates);

        for (unsigned int i=0; i<target_fractions.size(); ++i)
          Assert(cutoff_strain_rates[i] != std::numeric_limits<double>::max(),
                 ExcInternalError());

        constexpr unsigned int deformation_above_offset = 0;
        constexpr unsigned int area_above_offset = 2;
        constexpr unsigned int area_at_cutoff_offset = 4;

        double local_sums[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (const auto &entry : local_deformation_area_pairs)
          for (unsigned int i=0; i<target_fractions.size(); ++i)
            if (entry.first > cutoff_strain_rates[i])
              {
                local_sums[deformation_above_offset+i] += entry.first * entry.second;
                local_sums[area_above_offset+i] += entry.second;
              }
            else if (entry.first == cutoff_strain_rates[i])
              local_sums[area_at_cutoff_offset+i] += entry.second;

        double global_sums[6];
        Utilities::MPI::sum(local_sums, mpi_communicator, global_sums);

        std::array<double,2> area_fractions;
        for (unsigned int i=0; i<target_fractions.size(); ++i)
          {
            const double remaining_deformation =
              target_fractions[i] * global_total_deformation - global_sums[deformation_above_offset+i];
            const double deformation_at_cutoff =
              cutoff_strain_rates[i] * global_sums[area_at_cutoff_offset+i];
            const double tolerance = 100.0 * std::numeric_limits<double>::epsilon() * global_total_deformation;

            Assert(remaining_deformation >= -tolerance &&
                   remaining_deformation <= deformation_at_cutoff + tolerance,
                   ExcInternalError());

            const double area_at_cutoff =
              std::clamp(remaining_deformation / cutoff_strain_rates[i],
                         0.0,
                         global_sums[area_at_cutoff_offset+i]);

            area_fractions[i] =
              (global_sums[area_above_offset+i] + area_at_cutoff) / global_total_area;
          }

        return area_fractions;
      }
    }


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

      const std::array<double,2> area_fractions =
        compute_deformation_area_fractions(local_deformation_area_pairs,
                                           global_total_deformation,
                                           global_total_area,
                                           this->get_mpi_communicator());
      const double f80 = area_fractions[0];
      const double f90 = area_fractions[1];

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
