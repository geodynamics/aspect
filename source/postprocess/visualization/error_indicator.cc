// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2023 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/error_indicator.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      ErrorIndicator<dim>::
      ErrorIndicator ()
        :
        CellDataVectorCreator<dim>("")
      {}



      template <int dim>
      std::pair<std::string, std::unique_ptr<Vector<float>>>
      ErrorIndicator<dim>::execute() const
      {
        std::pair<std::string, std::unique_ptr<Vector<float>>>
        return_value ("error_indicator",
                      std::make_unique<Vector<float>>(this->get_triangulation().n_active_cells()));
        this->get_refinement_criteria(*return_value.second);

        return return_value;
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(ErrorIndicator,
                                                  "error indicator",
                                                  "A visualization output object that generates output "
                                                  "showing the estimated error or other mesh refinement "
                                                  "indicator as a spatially variable function with one "
                                                  "value per cell."
                                                  "\n\n"
                                                  "Physical units: None. (Strictly speaking, errors have "
                                                  "physical units of course, but because error "
                                                  "\\textit{indicators} can be computed from different "
                                                  "solution components and other input, we consider "
                                                  "error indicators unitless.)")

    }
  }
}
