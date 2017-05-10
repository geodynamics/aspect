/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/seismic_anomalies.h>
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
      std::pair<std::string, Vector<float> *>
      SeismicVsAnomaly<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("Vs_anomaly",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        // Calculate the Vs averages within n slices
        std::vector<double> Vs_depth_average(n_slices);
        this->get_lateral_averaging().get_Vs_averages(Vs_depth_average);

        // Estimates of the lateral Vs average at each depth are required
        // for all cell depths, including those
        // shallower than the midpoint of the first slice or
        // deeper than the midpoint of the last slice.
        // For this reason, the vector of Vs averages need to be padded by
        // one element before the first element and after the last element.
        // Mirror padding is chosen as it provides the least biased estimator
        // (in the absence of further information).
        std::vector<double> padded_Vs_depth_average(n_slices+2);
        padded_Vs_depth_average[0] = 2.*Vs_depth_average[0] - Vs_depth_average[1];
        padded_Vs_depth_average[n_slices+1] = 2.*Vs_depth_average[n_slices-1] - Vs_depth_average[n_slices-2];
        std::copy ( Vs_depth_average.begin(), Vs_depth_average.end(), padded_Vs_depth_average.begin() + 1 );



        // Calculate the maximum depth of the domain
        const double max_depth = this->get_geometry_model().maximal_depth();

        // The following lines evaluate the Vs at a single point per cell
        // and then compute the percentage deviation from the average in
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

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        // Loop over the cells
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              // Get the pressure, temperature and composition in the cell
              fe_values.reinit (cell);

              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                           in.temperature);
              fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                        in.pressure);
              fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                          in.velocity);
              fe_values[this->introspection().extractors.pressure].get_function_gradients (this->get_solution(),
                                                                                           in.pressure_gradient);
              in.position = fe_values.get_quadrature_points();

              // we do not need the strain rate
              in.strain_rate.resize(0);

              // Loop over compositional fields to get composition values
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    composition_values[c]);
              for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                {
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    in.composition[i][c] = composition_values[c][i];
                }
              in.cell = &cell;

              out.additional_outputs.push_back(
                std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
                (new MaterialModel::SeismicAdditionalOutputs<dim> (n_q_points)));
              this->get_material_model().evaluate(in, out);

              MaterialModel::SeismicAdditionalOutputs<dim> *seismic_outputs
                = out.template get_additional_output<MaterialModel::SeismicAdditionalOutputs<dim> >();
              const double Vs = seismic_outputs->vs[0];

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
              const double Vs_average = (1. - fractional_slice)*padded_Vs_depth_average[idx] + fractional_slice*padded_Vs_depth_average[idx+1];
              (*return_value.second)(cell_index) = (Vs - Vs_average)/Vs_average*1e2;
            }

        return return_value;
      }

      template <int dim>
      std::pair<std::string, Vector<float> *>
      SeismicVpAnomaly<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("Vp_anomaly",
                      new Vector<float>(this->get_triangulation().n_active_cells()));

        // Calculate the Vp averages within n slices
        std::vector<double> Vp_depth_average(n_slices);
        this->get_lateral_averaging().get_Vp_averages(Vp_depth_average);

        // Estimates of the lateral Vp average at each depth are required
        // for all cell depths, including those
        // shallower than the midpoint of the first slice or
        // deeper than the midpoint of the last slice.
        // For this reason, the vector of Vp averages need to be padded by
        // one element before the first element and after the last element.
        // Mirror padding is chosen as it provides the least biased estimator
        // (in the absence of further information).
        std::vector<double> padded_Vp_depth_average(n_slices+2);
        padded_Vp_depth_average[0] = 2.*Vp_depth_average[0] - Vp_depth_average[1];
        padded_Vp_depth_average[n_slices+1] = 2.*Vp_depth_average[n_slices-1] - Vp_depth_average[n_slices-2];
        std::copy ( Vp_depth_average.begin(), Vp_depth_average.end(), padded_Vp_depth_average.begin() + 1 );



        // Calculate the maximum depth of the domain
        const double max_depth = this->get_geometry_model().maximal_depth();

        // The following lines evaluate the Vp at a single point per cell
        // and then compute the percentage deviation from the average in
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

        std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

        // Loop over the cells
        typename DoFHandler<dim>::active_cell_iterator
        cell = this->get_dof_handler().begin_active(),
        endc = this->get_dof_handler().end();

        unsigned int cell_index = 0;
        for (; cell!=endc; ++cell,++cell_index)
          if (cell->is_locally_owned())
            {
              // Get the pressure, temperature and composition in the cell
              fe_values.reinit (cell);

              // get the various components of the solution, then
              // evaluate the material properties there
              fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                           in.temperature);
              fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(),
                                                                                        in.pressure);
              fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
                                                                                          in.velocity);
              fe_values[this->introspection().extractors.pressure].get_function_gradients (this->get_solution(),
                                                                                           in.pressure_gradient);
              in.position = fe_values.get_quadrature_points();

              // we do not need the strain rate
              in.strain_rate.resize(0);

              // Loop over compositional fields to get composition values
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values(this->get_solution(),
                    composition_values[c]);
              for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
                {
                  for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                    in.composition[i][c] = composition_values[c][i];
                }
              in.cell = &cell;

              out.additional_outputs.push_back(
                std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
                (new MaterialModel::SeismicAdditionalOutputs<dim> (n_q_points)));
              this->get_material_model().evaluate(in, out);

              MaterialModel::SeismicAdditionalOutputs<dim> *seismic_outputs
                = out.template get_additional_output<MaterialModel::SeismicAdditionalOutputs<dim> >();
              const double Vp = seismic_outputs->vp[0];

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
              const double Vp_average = (1. - fractional_slice)*padded_Vp_depth_average[idx] + fractional_slice*padded_Vp_depth_average[idx+1];
              (*return_value.second)(cell_index) = (Vp - Vp_average)/Vp_average*1e2;
            }

        return return_value;
      }

      template <int dim>
      void
      SeismicVsAnomaly<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Vs anomaly");
            {
              prm.declare_entry ("Number of depth slices", "50",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
                                 "average seismic shear wave velocities "
                                 "from which anomalies are calculated. "
                                 "Units: non-dimensional.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      SeismicVsAnomaly<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Vs anomaly");
            {
              n_slices = prm.get_integer ("Number of depth slices");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      SeismicVpAnomaly<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Vp anomaly");
            {
              prm.declare_entry ("Number of depth slices", "50",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
                                 "average seismic compressional wave velocities "
                                 "from which anomalies are calculated. "
                                 "Units: non-dimensional.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      SeismicVpAnomaly<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Vp anomaly");
            {
              n_slices = prm.get_integer ("Number of depth slices");
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SeismicVsAnomaly,
                                                  "Vs anomaly",
                                                  "A visualization output object that generates output "
                                                  "showing the percentage anomaly in the seismic "
                                                  "compressional wave speed $V_s$ as a spatially variable "
                                                  "function with one value per cell. This anomaly is shown "
                                                  "as a percentage change relative to the laterally averaged "
                                                  "velocity at the depth of the cell. This velocity is "
                                                  "calculated by linear interpolation between average values "
                                                  "calculated within equally thick depth slices. The "
                                                  "number of depth slices in the domain is user-defined. "
                                                  "Typically, the best results will be obtained if the number "
                                                  "of depth slices is balanced between being large enough to "
                                                  "capture step changes in velocities, but small enough to "
                                                  "maintain a reasonable number of evaluation points per slice. "
                                                  "Bear in mind that lateral averaging subsamples the "
                                                  "finite element mesh.")

      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SeismicVpAnomaly,
                                                  "Vp anomaly",
                                                  "A visualization output object that generates output "
                                                  "showing the percentage anomaly in the seismic "
                                                  "compressional wave speed $V_p$ as a spatially variable "
                                                  "function with one value per cell. This anomaly is shown "
                                                  "as a percentage change relative to the laterally averaged "
                                                  "velocity at the depth of the cell. This velocity is "
                                                  "calculated by linear interpolation between average values "
                                                  "calculated within equally thick depth slices. The "
                                                  "number of depth slices in the domain is user-defined. "
                                                  "Typically, the best results will be obtained if the number "
                                                  "of depth slices is balanced between being large enough to "
                                                  "capture step changes in velocities, but small enough to "
                                                  "maintain a reasonable number of evaluation points per slice. "
                                                  "Bear in mind that lateral averaging subsamples the "
                                                  "finite element mesh.")
    }
  }
}
