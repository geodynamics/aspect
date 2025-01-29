/*
  Copyright (C) 2014 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_prescribed_viscosity_h
#define _aspect_material_model_prescribed_viscosity_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parsed_function.h>


namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that applies a viscosity to a ''base model'' chosen from any of
     * the other available material models. This prescribed viscosity material model
     * allows the user to specify a function which describes where the viscosity should be
     * prescribed and a second function which describes the viscosity in that region.
     * This material model requires a base model which prescribes the viscosity and the
     * other material parameters in the rest of the model.
     * @ingroup MaterialModels
     */
    template <int dim>
    class PrescribedViscosity : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the base model at the beginning of the run.
         */
        void initialize() override;

        /**
         * Update the base model and viscosity function at the beginning of
         * each timestep.
         */
        void update() override;

        /**
         * Compute the material properties by evaluating the base model and
         * then overwrite the viscosity according to the viscosity function
         * in the locations set by the indicator function.
         */
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * Method to declare parameters related to prescribed viscosity model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to prescribed viscosity model
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Method that indicates whether material is compressible. Prescribed viscosity model is compressible
         * if and only if base model is compressible.
         */
        bool is_compressible () const override;

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        /**
         * Parsed function that specifies where the viscosity is prescribed.
         */
        Functions::ParsedFunction<dim> prescribed_viscosity_indicator_function;

        /**
         * Parsed function that specifies what the viscosity is set to in the
         * prescribed regions.
         */
        Functions::ParsedFunction<dim> prescribed_viscosity_function;

        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim>> base_model;
    };
  }
}

#endif
