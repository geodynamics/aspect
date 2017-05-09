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


#ifndef __aspect__newton_h
#define __aspect__newton_h

#include <aspect/simulator_access.h>
#include <aspect/global.h>
#include <aspect/assembly.h>
#include <aspect/material_model/interface.h>

#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  using namespace dealii;

  namespace MaterialModel
  {

    template <int dim>
    class MaterialModelDerivatives : public AdditionalMaterialOutputs<dim>
    {
      public:
        /**
         * Constructor. Initialize the various arrays of this structure with the
         * given number of quadrature points and (finite element) components.
         *
         * @param n_points The number of quadrature points for which input
         * quantities will be provided.
         * @param n_comp The number of vector quantities (in the order in which
         * the Introspection class reports them) for which input will be
         * provided.
         */
        MaterialModelDerivatives (const unsigned int n_points,
                                  const unsigned int /*n_comp*/)
        {
          dviscosities_dvelocity.resize(n_points, numbers::signaling_nan<double>());
          dviscosities_dpressure.resize(n_points, numbers::signaling_nan<double>());
          dviscosities_dtemperature.resize(n_points, numbers::signaling_nan<double>());
          dviscosities_dstrain_rate.resize(n_points, numbers::signaling_nan<SymmetricTensor<2,dim> >());
        };

        /**
         * The derivatives of the viscosities
         */
        std::vector<double> dviscosities_dvelocity;
        std::vector<double> dviscosities_dpressure;
        std::vector<double> dviscosities_dtemperature;
        std::vector<SymmetricTensor<2,dim> > dviscosities_dstrain_rate;
        std::vector<std::vector<double> > dviscosities_dcompositions;

    };
  }

  /**
   * Class that contains all runtime parameters and other helper functions
   * related to Newton solver. A global instance can be retrieved with
   * SimulatorAccess<dim>::get_newton_handler(), but keep in mind that it only
   * exists if parameters.include_melt_transport is true. TODO: rewrite
   */
  template <int dim>
  class NewtonHandler: public SimulatorAccess<dim>
  {
    public:
      NewtonHandler(ParameterHandler &prm);

      /**
       * Declare additional parameters that are needed in models with
       * melt transport (including the fluid pressure boundary conditions).
       */
      static void declare_parameters (ParameterHandler &prm);

      /**
       * Parse additional parameters that are needed in models with
       * melt transport (including the fluid pressure boundary conditions).
       *
       * This has to be called before edit_finite_element_variables,
       * so that the finite elements that are used for the additional melt
       * variables can be specified in the input file and are parsed before
       * the introspection object is created.
       */
      void parse_parameters (ParameterHandler &prm);

      /**
       * Create an additional material model output object that contains
       * the additional output variables needed in simulation with melt transport,
       * and attaches a pointer to it to the corresponding vector in the
       * MaterialModel::MaterialModelOutputs structure.
       */
      static void create_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &output);


  };

}

#endif
