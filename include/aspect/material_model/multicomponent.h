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

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A material model which is intended for use with multiple compositional fields.
     * For each material parameter the user supplies a comma delimited list of length N+1,
     * N is the number of compositional fields.  The additional field corresponds to the value
     * for background mantle.  They should be ordered ``background, composition1, composition2...''
     *
     * If a single value is given, then all the compositional fields are given that value. 
     * Other lengths of lists are not allowed.  For a given compositional field the material 
     * parameters are treated as constant, except density, which varies linearly with temperature
     * according to the thermal expansivity.. 
     *
     * When more than one field is present at a point, they are averaged arithmetically.
     * An exception is viscosity, which may be averaged arithmetically, harmonically, geometrically,
     * or by selecting the viscosity of the composition with the greatest volume fraction.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Multicomponent : public MaterialModel::InterfaceCompatibility<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * returns false
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
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

      private:
        void compute_volume_fractions( const std::vector<double> &compositional_fields, 
                                             std::vector<double> &fractions) const;
 
        double reference_T;

        enum {
            Harmonic,
            Arithmetic,
            Geometric,
            MaximumComposition
        } ViscosityAveraging;

        std::vector<double> densities;
        std::vector<double> viscosities;
        std::vector<double> thermal_expansivities;
        std::vector<double> thermal_conductivities;
        std::vector<double> specific_heats;
    };

  }
}

#endif
