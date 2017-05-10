/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#ifndef __aspect__material_model_force_vector_h
#define __aspect__material_model_force_vector_h

#include <aspect/material_model/interface.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional material model output that supplies a force vector.
     */
    template<int dim>
    class AdditionalMaterialOutputsForceVector: public AdditionalMaterialOutputs<dim>
    {
      public:
        AdditionalMaterialOutputsForceVector(const unsigned int n_points)
          : force_u(n_points), force_p(n_points), force_pc(n_points)
        {}

        virtual ~AdditionalMaterialOutputsForceVector()
        {}

        virtual void average (const MaterialAveraging::AveragingOperation /*operation*/,
                              const FullMatrix<double>  &/*projection_matrix*/,
                              const FullMatrix<double>  &/*expansion_matrix*/)
        {
          // TODO: not implemented
        }

        /**
         * Force tensor for the conservation of momentum equation (first
         * Stokes equation).
         */
        std::vector<Tensor<1,dim> > force_u;
        /**
         * Force for the conservation of mass equation (second Stokes
         * equation).
         */
        std::vector<double> force_p;
        /**
         * Force for the compaction pressure equation (when using melt
         * transport).
         */
        std::vector<double> force_pc;
    };

  }
}


#endif
