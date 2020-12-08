/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/postprocess/visualization/geoid.h>
#include <aspect/citation_info.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Geoid<dim>::
      Geoid ()
        :
        DataPostprocessorScalar<dim> ("geoid",
                                      update_quadrature_points)
      {}

      template <int dim>
      void
      Geoid<dim>::
      initialize()
      {
        CitationInfo::add("geoid");
      }

      template <int dim>
      void
      Geoid<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim> >(this->get_geometry_model()),
                     ExcMessage("The geoid postprocessor is currently only implemented for "
                                "the spherical shell geometry model."));

        for (unsigned int q=0; q<computed_quantities.size(); ++q)
          computed_quantities[q](0) = 0;

        const Postprocess::Geoid<dim> &geoid =
          this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::Geoid<dim> >();

#if DEAL_II_VERSION_GTE(9,3,0)
        auto cell = input_data.template get_cell<dim>();
#else
        auto cell = input_data.template get_cell<DoFHandler<dim> >();
#endif

        bool cell_at_top_boundary = false;
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->at_boundary(f) &&
              this->get_geometry_model().translate_id_to_symbol_name (cell->face(f)->boundary_id()) == "top")
            cell_at_top_boundary = true;

        if (cell_at_top_boundary)
          for (unsigned int q=0; q<input_data.evaluation_points.size(); ++q)
            computed_quantities[q](0) = geoid.evaluate(input_data.evaluation_points[q]);
      }

      template <int dim>
      std::list<std::string>
      Geoid<dim>::required_other_postprocessors() const
      {
        return std::list<std::string> (1, "geoid");
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Geoid,
                                                  "geoid",
                                                  "Visualization for the geoid solution. The geoid is given "
                                                  "by the equivalent water column height due to a gravity perturbation. "
                                                  "Units: \\si{\\meter}.")
    }
  }
}
