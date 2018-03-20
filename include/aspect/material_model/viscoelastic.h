/*
  Copyright (C) 2014 - 2017 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_viscoelastic_h
#define _aspect_material_model_viscoelastic_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model which is intended for simple linear viscoelastic
     * problems with multiple compositional fields.
     *
     * The viscoelastic constitutive relationship and implementation follows
     * a method commonly used in the Geodynamics community, where only the
     * elastic shear strength of rocks (e.g., shear modulus) is considered.
     * The material model directly follows the work flow and terminology
     * of Moresi et al. 2003 (J. Comp. Phys., v. 184, p. 476-497) equations
     * 23-32, which is commonly referred to in proceeding geodynamic
     * publications. However, a notable difference between this material
     * model and that of previous work is the use of compositional fields,
     * rather than tracers, to track and advect the stress tensor. The
     * material model will be updated when an option to track and advect
     * stresses with tracers is implemented.
     *
     * The first 3 (2D) or 6 (3D) compositional fields are used to track
     * individual components of the viscoelastic stress tensor, while additional
     * fields represent distinct rock types.
     *
     * The value of each compositional field representing distinct
     * rock types at a point is interpreted to be a volume fraction of that
     * rock type. If the sum of the compositional field volume fractions is
     * less than one, then the remainder of the volume is assumed to be
     * 'background material'.
     *
     * Several model parameters (densities, elastic shear moduli,
     * thermal expansivities, thermal conductivies, specific heats) can
     * be defined per-compositional field. For each material parameter the
     * user supplies a comma delimited list of length N+1, where N is the
     * number of compositional fields. The additional field corresponds to
     * the value for background material. They should be ordered
     * ``background, composition1, composition2...''. However, the first 3 (2D)
     * or 6 (3D) composition fields correspond to components of the elastic
     * stress tensor and their material values will not contribute to the volume
     * fractions. If a single value is given, then all the compositional fields
     * are given that value. Other lengths of lists are not allowed. For a given
     * compositional field the material parameters are treated as constant,
     * except density, which varies linearly with temperature according to the
     * thermal expansivity.
     *
     * When more than one compositional field is present at a point, they are
     * averaged arithmetically. An exception is viscosity, which may be averaged
     * arithmetically, harmonically, geometrically, or by selecting the
     * viscosity of the composition with the greatest volume fraction.
     *
     *
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Viscoelastic : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * The first 3 (2D) or 6 (3D) compositional fields are assumed
         * to be components of the viscoelastic stress tensor and
         * assigned volume fractions of zero.
         */
        const std::vector<double> compute_volume_fractions(
          const std::vector<double> &compositional_fields) const;

        /**
         * Reference temperature for thermal expansion.  All components use
         * the same reference_T.
         */
        double reference_T;

        /**
         * Enumeration for selecting which averaging scheme to use. Select
         * between harmonic, arithmetic, geometric, and maximum_composition.
         * The max composition scheme simply uses the parameter of whichever
         * field has the highest volume fraction.
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
         * Vector for field densities, read from parameter file.
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

        /**
         * Vector for field elastic shear moduli, read from parameter file.
         */
        std::vector<double> elastic_shear_moduli;

        /**
         * Bool indicating whether to use a fixed material time scale in the
         * viscoelastic rheology for all time steps (if true) or to use the
         * actual (variable) advection time step of the model (if false). Read
         * from parameter file.
         */
        bool use_fixed_elastic_time_step;

        /**
         * Bool indicating whether to use a stress averaging scheme to account
         * for differences between the numerical and fixed elastic time step
         * (if true). When set to false, the viscoelastic stresses are not
         * modified to account for differences between the viscoelastic time
         * step and the numerical time step. Read from parameter file.
         */
        bool use_stress_averaging;

        /**
         * Double for fixed elastic time step value, read from parameter file
         */
        double fixed_elastic_time_step;

    };

  }
}

#endif
