// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2024-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/depth_including_mesh_deformation.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/postprocess/current_surface.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DepthIncludingMeshDeformation<dim>::
      DepthIncludingMeshDeformation ()
        :
        DataPostprocessorScalar<dim> ("depth_including_mesh_deformation",
                                      update_quadrature_points),
        Interface<dim>("m")
      {}



      template <int dim>
      void
      DepthIncludingMeshDeformation<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();

        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        const Postprocess::CurrentSurface<dim> &surface =
          this->get_postprocess_manager().template get_matching_active_plugin<Postprocess::CurrentSurface<dim>>();

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = surface.depth_including_mesh_deformation (input_data.evaluation_points[q]);
      }



      template <int dim>
      std::list<std::string>
      DepthIncludingMeshDeformation<dim>::required_other_postprocessors() const
      {
        return {"current surface"};
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DepthIncludingMeshDeformation,
                                                  "depth including mesh deformation",
                                                  "A visualization output postprocessor that outputs "
                                                  "the depth for all points inside the domain, as "
                                                  "determined by the current model surface. This plugin "
                                                  "will include changes to the surface from mesh deformation."
                                                  "\n\n"
                                                  "Physical units: \\si{\\meter}.")
    }
  }
}
