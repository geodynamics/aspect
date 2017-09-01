/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#include <aspect/postprocess/heat_flux_densities.h>
#include <aspect/utilities.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    HeatFluxDensities<dim>::execute (TableHandler &statistics)
    {
      const char *unit = (dim==2)? "W/m" : "W/m^2";

      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim-1> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_gradients      | update_values |
                                        update_normal_vectors |
                                        update_q_points       | update_JxW_values);

      std::vector<Tensor<1,dim> > temperature_gradients (quadrature_formula.size());
      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      std::map<types::boundary_id, double> local_boundary_fluxes;
      std::map<types::boundary_id, double> local_areas;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      MaterialModel::MaterialModelInputs<dim> in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      // For every surface face on which it makes sense to compute a
      // heat flux and that is owned by this processor,
      // integrate the normal heat flux given by the formula
      //   j =  - k * n . grad T
      // Then compute the area of each boundary that lives on this processor.
      // Finally, sum over the processors and compute the ratio between the
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
              {
                fe_face_values.reinit (cell, f);
                // Set use_strain_rates to false since we don't need viscosity
                in.reinit(fe_face_values, cell, this->introspection(), this->get_solution(), false);

                this->get_material_model().evaluate(in, out);

                // Get the temperature gradients from the solution.
                fe_face_values[this->introspection().extractors.temperature].get_function_gradients (this->get_solution(), temperature_gradients);

                double local_normal_flux = 0;
                double local_area = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    const double thermal_conductivity
                      = out.thermal_conductivities[q];

                    local_normal_flux
                    +=
                      -thermal_conductivity *
                      (temperature_gradients[q] *
                       fe_face_values.normal_vector(q)) *
                      fe_face_values.JxW(q);

                    local_area += fe_face_values.JxW(q);
                  }

                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();
                local_boundary_fluxes[boundary_indicator] += local_normal_flux;
                local_areas[boundary_indicator] += local_area;
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_flux_densities;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();
        std::vector<double> local_boundary_fluxes_vector;
        std::vector<double> local_areas_vector;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p)
          {
            local_boundary_fluxes_vector.push_back (local_boundary_fluxes[*p]);
            local_areas_vector.push_back (local_areas[*p]);
          }

        // then collect contributions from all processors
        std::vector<double> global_values (local_boundary_fluxes_vector.size());
        std::vector<double> global_areas (local_areas_vector.size());
        Utilities::MPI::sum (local_boundary_fluxes_vector, this->get_mpi_communicator(), global_values);
        Utilities::MPI::sum (local_areas_vector, this->get_mpi_communicator(), global_areas);

        // and now take them apart into the global map again and compute ratios
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_flux_densities[*p] = global_values[index]/global_areas[index];
      }

      // now add all of the computed heat fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_flux_densities.begin();
           p != global_boundary_flux_densities.end(); ++p, ++index)
        {
          const std::string name = "Outward heat flux density for boundary with id "
                                   + Utilities::int_to_string(p->first)
                                   + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                 .translate_id_to_symbol_name (p->first))
                                   + " (" + unit +")";
          statistics.add_value (name, p->second);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name, 8);
          statistics.set_scientific (name, true);

          // finally have something for the screen
          screen_text.precision(4);
          screen_text << p->second << " " << unit
                      << (index == global_boundary_flux_densities.size()-1 ? "" : ", ");
        }

      return std::pair<std::string, std::string> ("Heat flux densities for boundaries:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(HeatFluxDensities,
                                  "heat flux densities",
                                  "A postprocessor that computes some statistics about "
                                  "the (conductive) heat flux density for each boundary "
                                  "id. The heat flux density is computed in outward "
                                  "direction, i.e., from the domain to the outside, using the "
                                  "formula $\\frac{1}{|\\Gamma_i|} \\int_{\\Gamma_i} -k \\nabla T \\cdot \\mathbf n$ "
                                  "where $\\Gamma_i$ is the part of the boundary with indicator $i$, "
                                  "$k$ is the thermal conductivity as reported by the material model, "
                                  "$T$ is the temperature, and $\\mathbf n$ is the outward normal. "
                                  "Note that the quantity so computed does not include any energy "
                                  "transported across the boundary by material transport in cases "
                                  "where $\\mathbf u \\cdot \\mathbf n \\neq 0$."
                                  "\n\n"
                                  "Note that the ``heat flux'' postprocessor computes the same "
                                  "quantity as the one here, but not divided by the area of "
                                  "the surface. In other words, it computes the "
                                  "\\textit{total} heat flux through each boundary."
                                 )
  }
}
