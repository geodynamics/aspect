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

#ifndef _aspect_material_model_additive_h
#define _aspect_material_model_additive_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that applies an average of the quadrature points in a cell  to
     * a ''base model'' chosen from any of the other available material models.
     * @ingroup MaterialModels
     */

    template <int dim>
    class Additive : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        virtual
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const;
        /**
         * Method to declare parameters related to depth-dependent model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to depth-dependent model
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * Method that indicates whether material is compressible. Depth dependent model is compressible
         * if and only if base model is compressible.
         */
        virtual bool is_compressible () const;

        /**
         * Method to calculate reference viscosity for the depth-dependent model. The reference
         * viscosity is determined by evaluating the depth-dependent part of the viscosity at
         * the mean depth of the model.
         */

        virtual double reference_viscosity () const;

        /**
         * Method to calculate reference density for the depth-dependent model. Because the depth-
         * dependent model deos not modify density, the reference density is equivalent to the
         * base model's reference density.
         */
        //virtual double reference_density () const;


      private:
        /**
         * Pointer to the material model used as the base model
         */
        std_cxx11::shared_ptr<MaterialModel::Interface<dim> > base_model;
    };
  }
}

#endif
