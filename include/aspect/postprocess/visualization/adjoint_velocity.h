/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_gravity_h
#define _aspect_postprocess_visualization_gravity_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that outputs the gravity.
       */
      template <int dim>
      class AdjointVelocity
        : public DataPostprocessorVector<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          AdjointVelocity ();

          virtual
          void
          compute_derived_quantities_vector (const std::vector<Vector<double> >              &solution_values,
                                             const std::vector<std::vector<Tensor<1,dim> > > &solution_gradients,
                                             const std::vector<std::vector<Tensor<2,dim> > > &solution_hessians,
                                             const std::vector<Point<dim> >                  &normals,
                                             const std::vector<Point<dim> >                  &evaluation_points,
                                             std::vector<Vector<double> >                    &computed_quantities) const;
      };
    }
  }
}

#endif
