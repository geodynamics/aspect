/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_entropy_model_h
#define _aspect_material_model_entropy_model_h

#include <aspect/material_model/interface.h>

#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/rheology/ascii_depth_profile.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model that is designed to use pressure and entropy (rather
     * than pressure and temperature) as independent variables. It will look up
     * all material properties in a data table for a given pressure and
     * entropy, and will additionally provide an additional output object of
     * type PrescribedTemperatureOutputs filled with the temperature, and
     * an additional output object of type PrescribedFieldOutput filled with
     * the densities (necessary for the projected density approximation of
     * the Stokes equation).
     * @ingroup MaterialModels
     */
    template <int dim>
    class EntropyModel: public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        virtual
        void
        initialize ();

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        virtual
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                 MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

        /**
         * Creates additional output objects of
         * type PrescribedTemperatureOutputs filled with the temperature
         * (necessary for solving the entropy equation), and
         * an additional output object of type PrescribedFieldOutput filled with
         * the densities (necessary for the projected density approximation of
         * the Stokes equation). Also creates SeismicAdditionalOutputs for
         * postprocessing purposes.
         */
        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


      private:
        /**
         * Minimum/Maximum viscosity and lateral viscosity variations.
         */
        double lateral_viscosity_prefactor;
        double min_eta;
        double max_eta;

        /**
         * The value for thermal conductivity. This model only
         * implements a constant thermal conductivity for the whole domain.
         */
        double thermal_conductivity_value;

        /**
         * Information about the location of data files.
         */
        std::string data_directory;
        std::string material_file_name;

        /**
         * Pointer to the StructuredDataLookup object that holds the material data.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<2>> material_lookup;

        /**
         * Pointer to the rheology model used for depth-dependence from an
         * ascii file
         */
        std::unique_ptr<Rheology::AsciiDepthProfile<dim>> depth_dependent_rheology;
    };
  }
}

#endif
