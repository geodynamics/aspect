// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2014-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_material_model_multicomponent_h
#define _aspect_material_model_multicomponent_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>


namespace aspect
{
  namespace MaterialModel
  {
    /**
     * An incompressible material model which is intended for use with multiple
     * compositional fields. Each compositional field is meant to be a single
     * rock type, where the value of the field at a point is interpreted to be
     * a volume fraction of that rock type.  If the sum of the compositional field
     * volume fractions is less than one, then the remainder of the volume is
     * assumed to be ``background mantle''.  If the sum of the compositional
     * field volume fractions is greater than one, then they are renormalized
     * to sum to one and there is no background mantle.
     *
     * For each material parameter the user supplies a comma delimited list of
     * length N+1, where N is the number of compositional fields.  The
     * additional field corresponds to the value for background mantle.  They
     * should be ordered ``background, composition1, composition2...''
     *
     * If a single value is given, then all the compositional fields are given
     * that value. Other lengths of lists are not allowed.  For a given
     * compositional field the material parameters are treated as constant,
     * except density, which varies linearly with temperature according to the
     * thermal expansivity.
     *
     * When more than one field is present at a point, they are averaged
     * arithmetically. An exception is viscosity, which may be averaged
     * arithmetically, harmonically, geometrically, or by selecting the
     * viscosity of the composition with the greatest volume fraction.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Multicomponent : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * This model is not compressible, so this returns false.
         */
        bool is_compressible () const override;

        /**
         * @}
         */

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
         * Reference temperature for thermal expansion.  All components use
         * the same reference_T.
         *
         * This variable is read from the parameter file through a parameter called 'Reference temperature'.
         */
        double reference_T;

        /**
         * Enumeration for selecting which viscosity averaging scheme to use.
         * This variable is read from the parameter file through a parameter called 'Viscosity averaging scheme'.
         */
        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging;

        /**
         * Vector for field viscosities, read from parameter file.
         * This variable is read from the parameter file through a parameter called 'Viscosities'.
         */
        std::vector<double> viscosities;

        /**
         * Vector for field thermal conductivities, read from parameter file.
         *
         * This variable is read from the parameter file through a parameter called 'Thermal conductivities'.
         */
        std::vector<double> thermal_conductivities;

        EquationOfState::MulticomponentIncompressible<dim> equation_of_state;
    };

  }
}

#endif
