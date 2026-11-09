// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/postprocess/visualization/partition.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Partition<dim>::
      Partition ()
        :
        // we don't need to know about any of the solution values
        // in order to determine the partition number. thus, no
        // need to specify any update flags
        DataPostprocessorScalar<dim> ("partition",
                                      update_default),
        Interface<dim>("")  // no physical units
      {}



      template <int dim>
      void
      Partition<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &/*input_data*/,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        Assert (computed_quantities[0].size() == 1, ExcInternalError());

        for (auto &quantity : computed_quantities)
          {
            // simply get the partition number from the triangulation
            quantity(0) = this->get_triangulation().locally_owned_subdomain();
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Partition,
                                                  "partition",
                                                  "A visualization output object that generates output "
                                                  "for the parallel partition that every cell of the "
                                                  "mesh is associated with."
                                                  "\n\n"
                                                  "Physical units: None.")
    }
  }
}
