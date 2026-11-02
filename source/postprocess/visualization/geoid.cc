// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/postprocess/visualization/geoid.h>

#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>

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
                                      update_quadrature_points),
        Interface<dim>("m")
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
                            std::vector<Vector<double>> &computed_quantities) const
      {
        AssertThrow (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()),
                     ExcMessage("The geoid postprocessor is currently only implemented for "
                                "the spherical shell geometry model."));

        for (auto &quantity : computed_quantities)
          quantity(0) = 0;

        const Postprocess::Geoid<dim> &geoid =
          this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::Geoid<dim>>();

        auto cell = input_data.template get_cell<dim>();

        bool cell_at_top_boundary = false;
        for (const unsigned int f : cell->face_indices())
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
        return {"geoid"};
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
                                                  "\n\n"
                                                  "Physical units: $\\text{m}$.")
    }
  }
}
