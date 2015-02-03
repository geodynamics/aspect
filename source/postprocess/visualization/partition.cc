/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/visualization/partition.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_out.h>


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
                                      update_default)
      {}



      template <int dim>
      void
      Partition<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> > &,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> > &,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        Assert (computed_quantities[0].size() == 1, ExcInternalError());

        for (unsigned int q=0; q<computed_quantities.size(); ++q)
          {
            // simply get the partition number from the triangulation
            computed_quantities[q](0) = this->get_triangulation().locally_owned_subdomain();
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
                                                  "mesh is associated with.")
    }
  }
}
