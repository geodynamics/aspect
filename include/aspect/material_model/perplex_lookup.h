/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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
        std::string perplex_file_name;
        double eta;
        double k_value;
        double min_temperature;
        double max_temperature;
        double min_pressure;
        double max_pressure;
    };

  }
}

#endif
