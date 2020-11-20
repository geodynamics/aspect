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


#include <aspect/postprocess/heat_flux_statistics.h>
#include <aspect/postprocess/heat_flux_map.h>

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
    HeatFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      std::vector<std::vector<std::pair<double, double> > > heat_flux_and_area =
        internal::compute_heat_flux_through_boundary_faces (*this);

      std::map<types::boundary_id, double> local_boundary_fluxes;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))
              {
                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();
                local_boundary_fluxes[boundary_indicator] += heat_flux_and_area[cell->active_cell_index()][f].first;
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
          local_values.emplace_back (local_boundary_fluxes[p]);

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

      // now add all of the computed heat fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_fluxes.begin();
           p != global_boundary_fluxes.end(); ++p, ++index)
        {
          const std::string name = "Outward heat flux through boundary with indicator "
                                   + Utilities::int_to_string(p->first)
                                   + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                 .translate_id_to_symbol_name (p->first))
                                   + " (W)";
          statistics.add_value (name, p->second);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name, 8);
          statistics.set_scientific (name, true);

          // finally have something for the screen
          screen_text.precision(4);
          screen_text << p->second << " W"
                      << (index == global_boundary_fluxes.size()-1 ? "" : ", ");
        }

      return std::pair<std::string, std::string> ("Heat fluxes through boundary parts:",
                                                  screen_text.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(HeatFluxStatistics,
                                  "heat flux statistics",
                                  "A postprocessor that computes some "
                                  "statistics about the heat flux "
                                  "density across each boundary in outward "
                                  "direction, i.e., from the domain to the "
                                  "outside. The heat flux is computed as sum "
                                  "of advective heat flux and conductive heat "
                                  "flux through Neumann boundaries, both "
                                  "computed as integral over the boundary area, "
                                  "and conductive heat flux through Dirichlet "
                                  "boundaries, which is computed using the "
                                  "consistent boundary flux method as described "
                                  "in ``Gresho, Lee, Sani, Maslanik, Eaton (1987). "
                                  "The consistent Galerkin FEM for computing "
                                  "derived boundary quantities in thermal and or "
                                  "fluids problems. International Journal for "
                                  "Numerical Methods in Fluids, 7(4), 371-394.''"
                                  "The point-wise heat flux can be obtained from the heat flux map postprocessor, "
                                  "which outputs the heat flux to a file, or the heat flux map "
                                  "visualization postprocessor, which outputs the heat flux for "
                                  "visualization. "
                                  "\n\n"
                                  "As stated, this postprocessor computes the \\textit{outbound} heat "
                                  "flux. If you "
                                  "are interested in the opposite direction, for example from "
                                  "the core into the mantle when the domain describes the "
                                  "mantle, then you need to multiply the result by -1."
                                  "\n\n"
                                  "\\note{In geodynamics, the term ``heat flux'' is often understood "
                                  "to be the quantity $- k \\nabla T$, which is really a heat "
                                  "flux \\textit{density}, i.e., a vector-valued field. In contrast "
                                  "to this, the current postprocessor only computes the integrated "
                                  "flux over each part of the boundary. Consequently, the units of "
                                  "the quantity computed here are $W=\\frac{J}{s}$.}"
                                  "\n\n"
                                  "The ``heat flux densities'' postprocessor computes the same "
                                  "quantity as the one here, but divided by the area of "
                                  "the surface.")
  }
}
