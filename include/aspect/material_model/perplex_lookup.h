// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_perplex_lookup_h
#define _aspect_material_model_perplex_lookup_h

#include <aspect/material_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that calls the thermodynamic software PerpleX
     * in order to evaluate material properties at a given point, namely
     * the density, heat capacity, thermal expansivity and compressibility.
     * The viscosity and thermal conductivity are globally constant.
     *
     * The usual PerpleX input files are required in the working directory:
     * an endmember dataset file, solution model dataset file and
     * option file. The path to a PerpleX input file must be given by
     * the user in the ASPECT parameter file; this file may be in any
     * directory. PerpleX output files will be written to the same
     * directory.
     *
     * Compositional fields are used to define the bulk composition at
     * each point. These compositional fields correspond to the
     * components as given in the PerpleX input file. The order of
     * components is preserved.
     *
     * This model requires PerpleX. A script to download and setup the
     * required files can be found in contrib/perplex. If the default
     * installation directory is not changed, cmake will automatically
     * find the required files during creation of the ASPECT build files.
     * See ./contrib/perplex/README.md
     *
     * WARNING: This model is extremely slow because there are many
     * redundant calls to evaluate the material properties;
     * it serves only as a proof of concept.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class PerpleXLookup : public Interface<dim>
    {
      public:
        /**
         * Initialization function. Loads the material data and sets up
         * pointers.
         */
        void
        initialize () override;

        bool is_compressible () const override;

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;


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
        void
        parse_parameters (ParameterHandler &prm) override;
        /**
         * @}
         */


      private:
        /**
         * This variable is read from the parameter file through a parameter called 'PerpleX input file name'.
         */
        std::string perplex_file_name;
        /**
         * This variable is read from the parameter file through a parameter called 'Viscosity'.
         */
        double eta;
        /**
         * This variable is read from the parameter file through a parameter called 'Thermal conductivity'.
         */
        double k_value;
        /**
         * This variable is read from the parameter file through a parameter called 'Minimum material temperature'.
         */
        double min_temperature;
        /**
         * This variable is read from the parameter file through a parameter called 'Maximum material temperature'.
         */
        double max_temperature;
        /**
         * This variable is read from the parameter file through a parameter called 'Minimum material pressure'.
         */
        double min_pressure;
        /**
         * This variable is read from the parameter file through a parameter called 'Maximum material pressure'.
         */
        double max_pressure;
    };

  }
}

#endif
