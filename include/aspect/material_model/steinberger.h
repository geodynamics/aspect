/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__model_steinberger_h
#define __aspect__model_steinberger_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace internal
    {
      class MaterialLookup;
      class LateralViscosityLookup;
      class RadialViscosityLookup;
    }
    /**
     * A variable viscosity material model that reads the essential values of
     * coefficients from tables in input files.
     *
     * The viscosity of this model is based on the paper
     * Steinberger/Calderwood 2006: "Models of large-scale viscous flow in the
     * Earth's mantle with contraints from mineral physics and surface
     * observations". The thermal conductivity is constant and the other
     * parameters are provided via lookup tables from the software PERPLEX.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Steinberger: public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * Called at the beginning of each time step and allows the material
         * model to update internal data structures.
         */
        virtual void update();
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

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double seismic_Vp (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const;

        virtual double seismic_Vs (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return true if the viscosity() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the density() function returns something that may
         * depend on the variable identifies by the argument.
         */
        virtual bool
        density_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the compressibility() function returns something
         * that may depend on the variable identifies by the argument.
         *
         * This function must return false for all possible arguments if the
         * is_compressible() function returns false.
         */
        virtual bool
        compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the specific_heat() function returns something that
         * may depend on the variable identifies by the argument.
         */
        virtual bool
        specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return true if the thermal_conductivity() function returns
         * something that may depend on the variable identifies by the
         * argument.
         */
        virtual bool
        thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
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
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual
        void
        evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                 typename Interface<dim>::MaterialModelOutputs &out) const;

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
        bool interpolation;
        bool latent_heat;
        bool compressible;
        bool use_lateral_average_temperature;
        double reference_eta;
        std::vector<double> avg_temp;
        double min_eta;
        double max_eta;
        double max_lateral_eta_variation;
        std::string datadirectory;
        std::vector<std::string> material_file_names;
        unsigned int n_material_data;
        std::string radial_viscosity_file_name;
        std::string lateral_viscosity_file_name;

        /**
         * In the incompressible case we need to adjust the temperature as if
         * there would be an adiabatic temperature increase to look up the
         * material properties in the lookup table.
         */
        double get_corrected_temperature (const double temperature,
                                          const double pressure,
                                          const Point<dim> &position) const;

        /**
         * In the incompressible case we need to adjust the pressure as if
         * there would be an compressible adiabatic pressure increase to look
         * up the material properties in the lookup table. Unfortunately we do
         * not know the adiabatic pressure profile for the incompressible case
         * and therefore we do not know the dynamic pressure. The only
         * currently possible solution is to use the adiabatic pressure
         * profile only, neglecting dynamic pressure for material lookup in
         * this case. This is essentially similar to having a depth dependent
         * reference profile for all properties and modifying the profiles
         * only in temperature-dimension.
         */
        double get_corrected_pressure (const double temperature,
                                       const double pressure,
                                       const Point<dim> &position) const;

        /**
         * This function returns the compressible density derived from the
         * list of loaded lookup tables.
         */
        double get_compressible_density (const double temperature,
                                         const double pressure,
                                         const std::vector<double> &compositional_fields,
                                         const Point<dim> &position) const;

        /**
         * We need to correct the compressible density in the incompressible
         * case to an incompressible profile. This is done by dividing the
         * compressible density with the density at this pressure at adiabatic
         * temperature and multiplying with the surface adiabatic density.
         */
        double get_corrected_density (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        /**
         * List of pointers to objects that read and process data we get from
         * Perplex files. There is one pointer/object per compositional field
         * data provided.
         */
        std::vector<std_cxx1x::shared_ptr<internal::MaterialLookup> > material_lookup;

        /**
         * Pointer to an object that reads and processes data for the lateral
         * temperature dependency of viscosity.
         */
        std_cxx1x::shared_ptr<internal::LateralViscosityLookup> lateral_viscosity_lookup;

        /**
         * Pointer to an object that reads and processes data for the radial
         * viscosity profile.
         */
        std_cxx1x::shared_ptr<internal::RadialViscosityLookup> radial_viscosity_lookup;

    };
  }
}

#endif
