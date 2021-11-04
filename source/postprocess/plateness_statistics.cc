/*
  Copyright (C) 2011-2021 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <aspect/postprocess/plateness_statistics.h>
#include <aspect/global.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    PlatenessStatistics<dim>::execute (TableHandler &statistics)
    {
      const types::boundary_id top_boundary_id = this->get_geometry_model().translate_symbolic_boundary_name_to_id("top");

      // create a quadrature formula for the velocity.
      const QGauss<dim-1> quadrature_formula (this->introspection().polynomial_degree.velocities+1);
   
      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);
         
      std::vector<Tensor<1,dim>> velocities (fe_face_values.n_quadrature_points);

      // Vector to store the integral of the second invariant of the strain rate
      // over the surface.
      double local_second_invariant_of_strain_rate_integral = 0.0;
      double local_surface_area_integral = 0.0;

      std::vector<Tensor<1,dim>> local_second_invariant_of_strain_rate(fe_face_values.n_quadrature_points);
      std::vector<Tensor<1,dim>> local_surface_area(fe_face_values.n_quadrature_points); 
      std::vector<std::pair<double,double>> local_second_invariant_of_strain_rate_and_corresponding_area; 
      std::vector<SymmetricTensor<2,dim>> strain_rate (quadrature_formula.size());
 
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned() && cell->at_boundary())
          {
            unsigned int face_idx = numbers::invalid_unsigned_int;
            bool at_upper_surface = false;
            {
              for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
                {
                  if (cell->at_boundary(f) && cell->face(f)->boundary_id() == top_boundary_id)
                    {
                      // If the cell is at the top boundary, assign face_idx.
                      face_idx = f;
                      at_upper_surface = true;
                      break;
                    }
                }

            }
            if (at_upper_surface)
              {
                fe_face_values.reinit(cell, face_idx);
                fe_face_values[this->introspection().extractors.velocities].get_function_symmetric_gradients (this->get_solution(), strain_rate);
              
                for (unsigned int q = 0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    //Merge so column 1 is strain rate invariants and column 2 are the corresponding areas
                    local_second_invariant_of_strain_rate[q] = std::max(std::sqrt(std::fabs(second_invariant(deviator(strain_rate[q])))),min_strain_rate);
                    local_surface_area[q] = fe_face_values.JxW(q);
                    
                    // second invariant of strain rate over the whole surface
                    local_second_invariant_of_strain_rate_integral = local_second_invariant_of_strain_rate[q] * local_surface_area[q];
                    local_surface_area_integral += local_surface_area[q];
                  }
              }
          }

      const double total_second_invariant_of_strain_rate = Utilities::MPI::sum (local_second_invariant_of_strain_rate_integral, this->get_mpi_communicator());
      const double total_surface_area = Utilities::MPI::sum (local_surface_area_integral, this->get_mpi_communicator());

      int n = sizeof(local_second_invariant_of_strain_rate)/sizeof(local_second_invariant_of_strain_rate[0]);
      // Entering values in vector of pairs
      for (int i=0; i<n; i++)
        local_second_invariant_of_strain_rate_and_corresponding_area.push_back( std::make_pair(local_second_invariant_of_strain_rate[i],local_surface_area[i]));
      //sort local second invariants of strain rate in size order
      std::sort(local_second_invariant_of_strain_rate_and_corresponding_area.begin(), local_second_invariant_of_strain_rate_and_corresponding_area.end());

      //find the fraction of the surface area where 80% of the deformation occurs
      //double f_80 = 0.0;
      double cumulative_second_invariant_of_strain_rate = 0.0;
      double cumulative_surface_area = 0.0;
      int row_number = 0;
      do
        {
          Assert(row_number < local_second_invariant_of_strain_rate_and_area.size(), ExcInternalError());
          cumulative_second_invariant_of_strain_rate = cumulative_second_invariant_of_strain_rate + local_second_invariant_of_strain_rate_and_corresponding_area[row_number].first * local_second_invariant_of_strain_rate_and_corresponding_area[row_number].second;
          cumulative_surface_area = cumulative_surface_area + local_second_invariant_of_strain_rate_and_corresponding_area[row_number].second;
          row_number = row_number + 1;
        }
      while (cumulative_second_invariant_of_strain_rate<0.8*total_second_invariant_of_strain_rate);
      double f_80 = cumulative_surface_area/total_surface_area;

      // Compute plateness
      // As explained in Tackely (2000) and Lourenco et al., 2020 for an isoviscous, 
      // internally heated calculation with Ra~10^6, 80% of the surface deformation
      // occurs in 60% of the surface area. Therefore we express plateness as described below:
      double plateness = 0.0;
      plateness = 1 - f_80/0.6;

      statistics.add_value ("Plateness ",
                            plateness);

      const std::string column = "Plateness ";
      statistics.set_precision (column, 8);
      statistics.set_scientific (column, true);

      std::ostringstream output;
      output.precision(4);

      output << plateness;

      return std::pair<std::string, std::string> ("Plateness :",
                                                  output.str());
    }



    template <int dim>
    void
    PlatenessStatistics<dim>::parse_parameters (ParameterHandler &prm)
   {
     min_strain_rate = prm.get_double("Minimum strain rate");
   }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(PlatenessStatistics,
                                  "plateness statistics",
                                  "A postprocessor that outputs the plateness. "
                                  "We define plateness following lourenco et al., 2020 and "
                                  "based on the works by Weinstein and Olson (1992) "
                                  "and Tackley (2000). The square root of the second invariant of "
                                  "strain rate, $\\dot{\\sigma_surf} = \\frac{\\dot{\\sigma_\\phi\\phi}}{\\sqrt(2)}$ "
                                  "is used. The total integrated $\\dot{\\sigma_surf}$ is calculated, "
                                  "followed by the fraction of the surface area in which 80\\% of "
                                  "that deformation occurs. This area fraction, denoted $f_{80}$ "
                                  "would be 0 for perfect plates (all deformation takes place "
                                  "within infinitely narrow zones). Tackley (2000) found that in "
                                  "isoviscous, internally heated calculations with Ra=$10^6$, "
                                  "$f_{80}\\approx$0.6. Therefore plateness should vary between "
                                  "0 for homogenous-viscosity cases and 1 for perfect plates, "
                                  "where plateness = $1 - \\frac{f_{80}}{0.6}$.")
  }
}
