/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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
