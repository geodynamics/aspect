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
      namespace
      {
        /**
         * Computes a simple 2n+1 point moving average of a vector.
         * http://en.wikipedia.org/wiki/Moving_average
	 * 
	 * Overwrites a vector<double> called values 
	 * with the moving average.
         */
        void
        compute_running_average(std::vector<double> &values,
                                const int ncells)
        {
          std::vector<double> temp(values.size());

	  // This nested loop runs over each of the values in the
	  // outer loop, and then sums the values within
	  // ncells of the value of interest within the inner loop.
	  // At the edges, the inner loop pads the array using the
	  // values at the edge of the array.
          for (unsigned int idx=0; idx<values.size(); idx++)
            {
              double sum = 0;
              for (int isum=-ncells; isum<=ncells; isum++)
                {
		  //
                  sum += values[std::max(0,std::min((int) values.size()-1,(int) idx+isum))];
                }
              temp[idx] = sum/((double) (ncells*2+1));
            }
          values = temp;
        }

      }

      template <int dim>
      std::pair<std::string, Vector<float> *>
      SeismicVsAnomaly<dim>::execute() const
      {
        std::pair<std::string, Vector<float> *>
        return_value ("Vs_anomaly",
                      new Vector<float>(this->get_triangulation().n_active_cells()));
		
	// These three lines calculate the Vs averages within n slices,
	// and then calculates a simple (unweighted) running average
	// (effectively smoothing the averages)
	// The running average simply calculates the mean of the
	// slice of interest and the slices either side.
	
	std::vector<double> Vs_depth_average(n_slices);
        const unsigned int npoints = 2; 
        this->get_lateral_averaging().get_Vs_averages(Vs_depth_average);
        compute_running_average(Vs_depth_average, npoints);
	
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
	      // the depth slice within which that point resides
              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
              const unsigned int idx = static_cast<unsigned int>((depth*n_slices)/max_depth);
              Assert(idx<n_slices, ExcInternalError());

              // Compute the percentage deviation from the average
              (*return_value.second)(cell_index) = (Vs - Vs_depth_average[idx])/Vs_depth_average[idx]*1e2;
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
		
	// These three lines calculate the Vp averages within n slices,
	// and then calculates a simple (unweighted) running average
	// (effectively smoothing the averages)
	// The running average simply calculates the mean of the
	// slice of interest and the slices either side.
	std::vector<double> Vp_depth_average(n_slices);
        const unsigned int npoints = 2; 
        this->get_lateral_averaging().get_Vp_averages(Vp_depth_average);
        compute_running_average(Vp_depth_average, npoints);
	
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
	      // the depth slice within which that point resides
              const double depth = this->get_geometry_model().depth(fe_values.quadrature_point(0));
              const unsigned int idx = static_cast<unsigned int>((depth*n_slices)/max_depth);
              Assert(idx<n_slices, ExcInternalError());

              // Compute the percentage deviation from the average
              (*return_value.second)(cell_index) = (Vp - Vp_depth_average[idx])/Vp_depth_average[idx]*1e2;
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
            prm.enter_subsection("Seismic Vs anomaly");
            {
              prm.declare_entry ("Number of depth slices", "50",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
				 "average seismic shear wave velocities "
				 "from which anomalies are calculated. "
                                 "Units: non-dimensional.");
	    }
	  }
	}
      }
	      
      template <int dim>
      void
      SeismicVsAnomaly<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Seismic Vs anomaly");
            {
              n_slices = prm.get_integer ("Number of depth slices");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
	}
      }

      template <int dim>
      void
      SeismicVpAnomaly<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Seismic Vp anomaly");
            {
              prm.declare_entry ("Number of depth slices", "50",
                                 Patterns::Integer (1),
                                 "Number of depth slices used to define "
				 "average seismic compressional wave velocities "
				 "from which anomalies are calculated. "
                                 "Units: non-dimensional.");
	    }
	  }
	}
      }
	      
      template <int dim>
      void
      SeismicVpAnomaly<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Seismic Vp anomaly");
            {
              n_slices = prm.get_integer ("Number of depth slices");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
	}
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
						  "as a percentage change relative to the unweighted moving "
						  "average of laterally averaged velocities within the "
						  "three depth slices containing and adjacent to this cell.")

      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SeismicVpAnomaly,
                                                  "Vp anomaly",
                                                  "A visualization output object that generates output "
						  "showing the percentage anomaly in the seismic "
						  "compressional wave speed $V_p$ as a spatially variable "
						  "function with one value per cell. This anomaly is shown "
						  "as a percentage change relative to the unweighted moving "
						  "average of laterally averaged velocities within the "
						  "three depth slices containing and adjacent to this cell.")
    }
  }
}
