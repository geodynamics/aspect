/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

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
                                 const MPI_Comm &comm);

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
                                const MPI_Comm &comm);

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
        * Returns the volume fractions of all the phases over all of the current
        * evaluation points. The outer loop is over the phases, and the inner
        * loop is over all the evaluation points.
        */
        void fill_phase_volume_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                          NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out) const;

        /**
         * Returns the cell-wise averaged enthalpy derivatives for the evaluate
         * function and postprocessors. The function returns two pairs, the
         * first one represents the temperature derivative, the second one the
         * pressure derivative. The first member of each pair is the derivative,
         * the second one the number of vertex combinations the function could
         * use to compute the derivative. The second member is useful to handle
         * the case no suitable combination of vertices could be found (e.g.
         * if the temperature and pressure on all vertices of the current
         * cell is identical.
         */
        std::array<std::pair<double, unsigned int>,2>
        enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const;
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
         * @name Reference quantities
         * @{
         */
        double reference_viscosity () const override;
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
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
        bool interpolation;
        bool latent_heat;
        bool use_lateral_average_temperature;

        /**
         * Reference viscosity. Only used for pressure scaling purposes
         * and returned by the reference_viscosity() function.
         */
        double reference_eta;

        /**
         * The value for thermal conductivity. This model only
         * implements a constant thermal conductivity for the whole domain.
         */
        double thermal_conductivity_value;

        /**
         * Information about lateral temperature averages.
         */
        std::vector<double> avg_temp;
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
        std::vector<std::string> material_file_names;
        std::string radial_viscosity_file_name;
        std::string lateral_viscosity_file_name;

        /**
         * List of pointers to objects that read and process data we get from
         * Perplex files.
         */
        std::vector<std::unique_ptr<MaterialModel::MaterialUtilities::Lookup::PerplexReader> > material_lookup;

        /**
        * Vector of strings containing the names of the unique phases in all the material lookups.
        */
        std::vector<std::string> unique_phase_names;

        /**
        * Vector of vector of unsigned ints which constitutes mappings
        * between lookup phase name vectors and unique_phase_names.
        * The element unique_phase_indices[i][j] contains the
        * index of phase name j from lookup i as it is found in unique_phase_names.
        */
        std::vector<std::vector<unsigned int>> unique_phase_indices;

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

    };
  }
}

#endif
