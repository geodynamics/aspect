/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_steinberger_h
#define _aspect_material_model_steinberger_h

#include <aspect/material_model/interface.h>
#include <aspect/material_model/equation_of_state/thermodynamic_table_lookup.h>
#include <aspect/material_model/thermal_conductivity/interface.h>

#include <aspect/simulator_access.h>
#include <deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace internal
    {
      /**
       * A class that reads in a text file that contains the
       * temperature-dependency of viscosity for a set of equidistant depth
       * layers. See
       * the data/material-model/steinberger directory for an example data
       * file.
       * The class can return the value for a given depth.
       */
      class LateralViscosityLookup
      {
        public:
          /**
           * Read in a file.
           */
          LateralViscosityLookup(const std::string &filename,
                                 const MPI_Comm comm);

          /**
           * Returns a temperature-dependency for a given depth.
           */
          double lateral_viscosity(double depth) const;

          /**
           * Number of depth slices of the read file.
           */
          int get_nslices() const;
        private:
          /**
           * Stored values
           */
          std::vector<double> values;

          /**
           * Stored bounds an depths.
           */
          double min_depth;
          double delta_depth;
          double max_depth;
      };

      /**
       * A class that reads in a text file that contains the
       * viscosity for a set of equidistant depth layers. See
       * the data/material-model/steinberger directory for an example data
       * file.
       * The class can return the value for a given depth.
       */
      class RadialViscosityLookup
      {
        public:
          /**
           * Constructor. Reads in the given file.
           */
          RadialViscosityLookup(const std::string &filename,
                                const MPI_Comm comm);

          /**
           * Return the viscosity for a given depth.
           */
          double radial_viscosity(double depth) const;

        private:
          /**
           * Stored data values.
           */
          std::vector<double> values;

          /**
           * Depth bounds for the read in values.
           */
          double min_depth;
          double delta_depth;
          double max_depth;
      };
    }

    /**
     * A variable viscosity material model that reads the essential values of
     * coefficients from tables in input files.
     *
     * The viscosity of this model is based on the paper
     * Steinberger & Calderwood 2006: "Models of large-scale viscous flow in the
     * Earth's mantle with constraints from mineral physics and surface
     * observations". The thermal conductivity is constant or follows a
     * pressure-temperature dependent approximation and the other
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
        void
        initialize () override;

        /**
         * Called at the beginning of each time step and allows the material
         * model to update internal data structures.
         */
        void update() override;

        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in.
         */
        void
        evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
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

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;


      private:
        /**
         * Whether the compositional fields representing mass fractions
         * should be normalized to one when computing their fractions
         * (if false), or whether there is an additional composition
         * (the background field) that is not represented by a
         * compositional field, and makes up the remaining fraction of
         * material if the compositional fields add up to less than one
         * at any given location (if true).
         */
        bool has_background_field;

        /**
         * Pointer to a composition mask, which is meant to be filled with
         * one entry per compositional field that determines if this
         * field is considered to represent a mass fractions (if the entry
         * is set to true) or not (if set to false). This is needed for
         * averaging of material properties.
         */
        std::unique_ptr<ComponentMask> composition_mask;

        /**
         * The thermodynamic lookup equation of state.
         */
        EquationOfState::ThermodynamicTableLookup<dim> equation_of_state;

        /**
         * Boolean describing whether to use the lateral average temperature
         * for computing the viscosity, rather than the temperature
         * on the reference adiabat.
         */
        bool use_lateral_average_temperature;

        /**
         * The thermal conductivity parametrization to use. This material
         * model supports either a constant thermal conductivity or a
         * pressure- and temperature-dependent thermal conductivity.
         */
        std::unique_ptr<ThermalConductivity::Interface<dim>> thermal_conductivity;

        /**
         * Compositional prefactors with which to multiply the reference viscosity.
         * Volume fractions are used to weight the prefactors according to the
         * assigned viscosity averaging scheme.
         */
        std::vector<double> viscosity_prefactors;
        MaterialUtilities::CompositionalAveragingOperation viscosity_averaging_scheme;

        /**
         * Information about lateral temperature averages.
         */
        std::vector<double> average_temperature;
        unsigned int n_lateral_slices;

        /**
         * Minimum and maximum allowed viscosity, as well as the maximum allowed
         * viscosity variation compared to the average radial viscosity.
         */
        double min_eta;
        double max_eta;
        double max_lateral_eta_variation;

        /**
         * Information about the location of data files.
         */
        std::string data_directory;
        std::string radial_viscosity_file_name;
        std::string lateral_viscosity_file_name;

        /**
         * Pointer to an object that reads and processes data for the lateral
         * temperature dependency of viscosity.
         */
        std::unique_ptr<internal::LateralViscosityLookup> lateral_viscosity_lookup;

        /**
         * Pointer to an object that reads and processes data for the radial
         * viscosity profile.
         */
        std::unique_ptr<internal::RadialViscosityLookup> radial_viscosity_lookup;

        /**
         * A function that fills the prescribed additional outputs in the
         * MaterialModelOutputs object that is handed over, if it exists,
         * in this case, densities for the projected density approximation.
         * Does nothing otherwise.
         */
        void fill_prescribed_outputs (const unsigned int i,
                                      const std::vector<double> &volume_fractions,
                                      const MaterialModel::MaterialModelInputs<dim> &in,
                                      MaterialModel::MaterialModelOutputs<dim> &out) const;

    };
  }
}

#endif
