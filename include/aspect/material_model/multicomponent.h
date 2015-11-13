/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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


#ifndef __aspect__model_multicomponent_h
#define __aspect__model_multicomponent_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model which is intended for use with multiple compositional
     * fields. Each compositional field is meant to be a single rock type,
     * where the value of the field at a point is interpreted to be a volume
     * fraction of that rock type.  If the sum of the compositional field
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
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * This model is not compressible, so this returns false.
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

        double reference_thermal_diffusivity () const;

        double reference_cp () const;
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
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        /**
         * From a list of compositional fields of length N, we come up with an
         * N+1 length list that which also includes the fraction of
         * ``background mantle''. This list should sum to one, and is
         * interpreted as volume fractions.  If the sum of the
         * compositional_fields is greater than one, we assume that there is
         * no background mantle (i.e., that field value is zero).  Otherwise,
         * the difference between the sum of the compositional fields and 1.0
         * is assumed to be the amount of background mantle.
         */
        const std::vector<double> compute_volume_fractions(
          const std::vector<double> &compositional_fields) const;
        /**
         * Reference temperature for thermal expansion.  All components use
         * the same reference_T.
         */
        double reference_T;

        /**
        * Enumeration for selecting which averaging scheme to use.
        * Select between harmonic, arithmetic, geometric, and
        * maximum_composition.  The max composition scheme simply uses the
        * parameter of whichever field has the highest volume fraction.
        */
        enum AveragingScheme
        {
          harmonic,
          arithmetic,
          geometric,
          maximum_composition
        };


        AveragingScheme viscosity_averaging;

        double average_value (const std::vector<double> &composition,
                              const std::vector<double> &parameter_values,
                              const enum AveragingScheme &average_type) const;


        /**
         * Vector for field densities, read from parameter file .
         */
        std::vector<double> densities;

        /**
         * Vector for field viscosities, read from parameter file.
         */
        std::vector<double> viscosities;

        /**
         * Vector for field thermal expnsivities, read from parameter file.
         */
        std::vector<double> thermal_expansivities;

        /**
         * Vector for field thermal conductivities, read from parameter file.
         */
        std::vector<double> thermal_conductivities;

        /**
         * Vector for field specific heats, read from parameter file.
         */
        std::vector<double> specific_heats;
    };

  }
}

#endif
