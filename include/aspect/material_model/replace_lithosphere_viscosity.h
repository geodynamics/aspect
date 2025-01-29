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

#ifndef _aspect_material_model_replace_lithosphere_viscosity_h
#define _aspect_material_model_replace_lithosphere_viscosity_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/initial_temperature/lithosphere_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that applies a given constant viscosity in the lithosphere.
     * Viscosity below this is taken from a ''base model'' chosen from any of the
     * other available material models. The 'replace lithosphere viscosity'
     * material model allows the user to specify the depth of the lithosphere-asthenosphere
     * boundary either as one value or as a file. All other properties are derived
     * from the base model.
     * @ingroup MaterialModels
     */
    template <int dim>
    class ReplaceLithosphereViscosity : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Initialize the base model at the beginning of the run.
         */
        void initialize() override;

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * Method to declare parameters related to replace lithosphere viscosity model
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Method to parse parameters related to replace lithosphere viscosity model
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Method that indicates whether material is compressible. Replace lithosphere viscosity model is compressible
         * if and only if base model is compressible.
         */
        bool is_compressible () const override;

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        /**
         * This parameter gives the viscosity set within the lithosphere.
         */
        double lithosphere_viscosity;


        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim>> base_model;


        InitialTemperature::LABDepth::LABDepthLookup<dim> lab_depth_lookup;
    };
  }
}

#endif
