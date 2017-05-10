/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#ifndef __aspect__newton_h
#define __aspect__newton_h

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/material_model/interface.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  using namespace dealii;

  namespace MaterialModel
  {
    /**
     * This class holds the derivatives for the Newton solver.
     */
    template <int dim>
    class MaterialModelDerivatives : public AdditionalMaterialOutputs<dim>
    {
      public:
        /**
         * Constructor. Initialize the various arrays of this structure with the
         * given number of quadrature points.
         */
        MaterialModelDerivatives (const unsigned int n_points);

        /**
         * The derivatives of the viscosities
         */
        std::vector<double> viscosity_derivative_wrt_pressure;
        std::vector<SymmetricTensor<2,dim> > viscosity_derivative_wrt_strain_rate;

    };
  }

  /**
   * A Class which can declare and parse parameters and creates
   * material model outputs for the Newton solver.
   */
  template <int dim>
  class NewtonHandler: public SimulatorAccess<dim>
  {
    public:
      /**
       * Create an additional material model output object that contains
       * the additional output variables (the derivatives) needed for the
       * Newton solver.
       */
      static void create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output);
  };

}

#endif
