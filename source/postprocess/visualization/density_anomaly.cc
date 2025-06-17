/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/density_anomaly.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/lateral_averaging.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DensityAnomaly<dim>::
      DensityAnomaly ()
        :
        CellDataVectorCreator<dim>("kg/m^3")
      {}



      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      DensityAnomaly<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>>
        return_value ("density_anomaly",
                      std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));

        std::vector<double> padded_density_depth_average;
        if (average_density_scheme == lateral_average)
          {
            std::vector<double> density_depth_average(n_slices);
            this->get_lateral_averaging().get_density_averages(density_depth_average);

            // Estimates of the lateral density average at each depth are required
            // for all cell depths, including those
            // shallower than the midpoint of the first slice or
            // deeper than the midpoint of the last slice.
            // This way, the depth_average_density used for each cell is interpolated
            // from the padded_density_depth_average right above and below the cell.
            //
            // ---cell(n-1)---cell(n)---cell(n+1)---
            //                   |
            //                   |
            //        (fractional|
            //            _slice)|
            // --slice(idx)-----------slice(idx+1)----
            padded_density_depth_average.resize(n_slices+2);

            padded_density_depth_average[0] = 2.*density_depth_average[0] - density_depth_average[1];
            padded_density_depth_average[n_slices+1] = 2.*density_depth_average[n_slices-1] - density_depth_average[n_slices-2];

            std::copy (density_depth_average.begin(), density_depth_average.end(), padded_density_depth_average.begin() + 1 );
          }

        // The following lines evaluate the density at a single point per cell
        // and then compute the difference from the average in
        // the slice within which the point lies
        const QMidpoint<dim> quadrature_formula;
        const unsigned int n_q_points = quadrature_formula.size();


        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values   |
                                 update_gradients |
                                 update_quadrature_points );

        MaterialModel::MaterialModelInputs<dim> in(n_q_points, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(n_q_points, this->n_compositional_fields());
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              // Get the pressure, temperature and composition in the cell
              fe_values.reinit (cell);

              in.reinit(fe_values, cell, this->introspection(), this->get_solution());

              // Evaluate density at the cell center
              this->get_material_model().evaluate(in, out);
              const double density = out.densities[0];

              // Compute the absolute deviation from the average
              switch (average_density_scheme)
                {
                  case reference_profile:
                  {
                    // Evaluate adiabatic/reference density at the cell center
                    const double adiabatic_density = this->get_adiabatic_conditions().density(in.position[0]);
                    (*return_value.second)(cell->active_cell_index()) = density - adiabatic_density;
                    break;
                  }
                  case lateral_average:
                  {
                    // Calculate the maximum depth of the domain
                    const double max_depth = this->get_geometry_model().maximal_depth();

                    // Find the depth of the zeroth quadrature point in the cell and work out
                    // the depth slice which has its center just above that point
                    const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));

                    // Please note: static_cast<int> always truncates the double (i.e. 1.9 -> 1)
                    // Here we subtract 0.5 to find the first slice with its center shallower than the cell,
                    // and add 1 to take padding into account (-0.5 + 1.0 = +0.5)
                    const double slice_depth = (depth*n_slices)/max_depth + 0.5;
                    const unsigned int idx = static_cast<unsigned int>(slice_depth);
                    const double fractional_slice = slice_depth - static_cast<double>(idx);
                    Assert(idx<n_slices+1, ExcInternalError());

                    // Compute the percentage deviation from the average
                    const double depth_average_density = (1. - fractional_slice)*padded_density_depth_average[idx] + fractional_slice*padded_density_depth_average[idx+1];
                    (*return_value.second)(cell->active_cell_index()) = density - depth_average_density;
                    break;
                  }
                }
            }
        return return_value;
      }



      template <int dim>
      void
      DensityAnomaly<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Density anomaly");
            {
              prm.declare_entry ("Average density scheme", "reference profile",
                                 Patterns::Selection("reference profile|lateral average"),
                                 "Scheme to compute the average density-depth "
                                 "profile. The reference profile option evaluates "
                                 "the conditions along the reference adiabat "
                                 "according to the material model. "
                                 "The lateral average option instead calculates "
                                 "a lateral average from subdivision of the mesh. "
                                 "The lateral average option may produce spurious "
                                 "results where there are sharp density changes.");

              prm.declare_entry ("Number of depth slices","20",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
                                 "average density.");

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      DensityAnomaly<dim>::parse_parameters(ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Density anomaly");
            {
              // Average density scheme
              if (prm.get ("Average density scheme") == "reference profile")
                average_density_scheme = reference_profile;
              else if (prm.get ("Average density scheme") == "lateral average")
                average_density_scheme = lateral_average;
              else
                AssertThrow(false, ExcMessage("Not a valid average density scheme."));

              n_slices = prm.get_integer("Number of depth slices");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DensityAnomaly,
                                                  "density anomaly",
                                                  "A visualization output postprocessor that outputs the density minus the depth-average of the density."
                                                  "In the ``lateral average'' scheme, the average density is calculated using the lateral averaging function"
                                                  "from the ``depth average'' postprocessor and interpolated linearly between the layers specified through "
                                                  "``Number of depth slices''. In the ``reference profile'' scheme, the adiabatic density is used as the"
                                                  "average density."
                                                  "\n\n"
                                                  "Physical units: \\si{\\kg/m^3}.")
    }
  }
}
