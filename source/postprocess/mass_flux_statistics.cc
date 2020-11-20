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


#include <aspect/postprocess/mass_flux_statistics.h>
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
    MassFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      // First determine the units for the output
      const std::string unit = (this->convert_output_to_years())
                               ?
                               "kg/yr"
                               :
                               "kg/s";
      const double in_years = (this->convert_output_to_years())
                              ?
                              year_in_seconds
                              :
                              1.0;

      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim-1> quadrature_formula (this->introspection().polynomial_degree.velocities + 1);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_values            | update_gradients |
                                        update_normal_vectors    |
                                        update_quadrature_points | update_JxW_values);

      std::vector<std::vector<double> > composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      std::map<types::boundary_id, double> local_boundary_fluxes;

      MaterialModel::MaterialModelInputs<dim> in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      // for every surface face on which it makes sense to compute a
      // mass flux and that is owned by this processor,
      // integrate the normal mass flux given by the formula
      //   j =  \rho * v * n
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
              {
                fe_face_values.reinit (cell, f);
                // Set use_strain_rates to false since we don't need viscosity
                in.reinit(fe_face_values, cell, this->introspection(), this->get_solution(), false);

                this->get_material_model().evaluate(in, out);


                double local_normal_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                    local_normal_flux
                    +=
                      out.densities[q]
                      * (in.velocity[q] * fe_face_values.normal_vector(q))
                      * fe_face_values.JxW(q);
                  }

                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();
                local_boundary_fluxes[boundary_indicator] += local_normal_flux * in_years;
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_fluxes;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();
        std::vector<double> local_values;
        local_values.reserve(boundary_indicators.size());
        for (const auto p : boundary_indicators)
          local_values.push_back (local_boundary_fluxes[p]);

        // then collect contributions from all processors
        std::vector<double> global_values (local_values.size());
        Utilities::MPI::sum (local_values, this->get_mpi_communicator(), global_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes[*p] = global_values[index];
      }

      // now add all of the computed mass fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_fluxes.begin();
           p != global_boundary_fluxes.end(); ++p, ++index)
        {
          const std::string name = "Outward mass flux through boundary with indicator "
                                   + Utilities::int_to_string(p->first)
                                   + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                 .translate_id_to_symbol_name (p->first))
                                   + " (" + unit + ")";
          statistics.add_value (name, p->second);

          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name, 8);
          statistics.set_scientific (name, true);

          // finally have something for the screen
          screen_text.precision(4);
          screen_text << p->second << " " << unit
                      << (index == global_boundary_fluxes.size()-1 ? "" : ", ");
        }

      return std::pair<std::string, std::string> ("Mass fluxes through boundary parts:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MassFluxStatistics,
                                  "mass flux statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the mass flux across boundaries. For each boundary "
                                  "indicator (see your geometry description for which boundary "
                                  "indicators are used), the mass flux is computed in outward "
                                  "direction, i.e., from the domain to the outside, using the "
                                  "formula $\\int_{\\Gamma_i} \\rho \\mathbf v \\cdot \\mathbf n$ "
                                  "where $\\Gamma_i$ is the part of the boundary with indicator $i$, "
                                  "$\\rho$ is the density as reported by the material model, "
                                  "$\\mathbf v$ is the velocity, and $\\mathbf n$ is the outward normal. "
                                  "\n\n"
                                  "As stated, this postprocessor computes the \\textit{outbound} mass "
                                  "flux. If you "
                                  "are interested in the opposite direction, for example from "
                                  "the core into the mantle when the domain describes the "
                                  "mantle, then you need to multiply the result by -1."
                                  "\n\n"
                                  "\\note{In geodynamics, the term ``mass flux'' is often understood "
                                  "to be the quantity $\\rho \\mathbf v$, which is really a mass "
                                  "flux \\textit{density}, i.e., a vector-valued field. In contrast "
                                  "to this, the current postprocessor only computes the integrated "
                                  "flux over each part of the boundary. Consequently, the units of "
                                  "the quantity computed here are $\\frac{kg}{s}$.}")
  }
}
