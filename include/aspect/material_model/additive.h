/*
  Copyright (C) 2021 by the authors of the ASPECT code.

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
     * A material model that uses the current compositional field to change the
     * the starting viscosity and density field of the ''base model''. This base
     * model can be chosen from any of the other available material models
     * @ingroup MaterialModels
     */
    template <int dim>
    class Additive : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the base model at the beginning of the run.
         */
        virtual
        void initialize();

        /**
         * Update the base model at the beginning of each timestep.
         */
        virtual
        void update();

        /**
         * Function to update the density and viscosity using the values
         * of the compositional field
         */
        virtual
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const;

        /**
         * Method to declare parameters related to additive model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to additive model
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * Method that indicates whether the material is compressible. The additive model is compressible
         * if and only if base model is compressible.
         */
        virtual bool is_compressible () const;

        /**
         * Method to calculate reference viscosity for the additive model. The reference
         * viscosity is determined by evaluating the viscosity at
         * the mean depth of the model.
         */
        virtual double reference_viscosity () const;


      private:
        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim> > base_model;

        /**
        * Value with which the density update will be multiplied in each adjoint iteration
         */
        double density_update_factor;

        /**
         * Value with which the viscosity update will be multiplied in each adjoint iteration
         */
        double viscosity_update_factor;

    };
  }
}

#endif
