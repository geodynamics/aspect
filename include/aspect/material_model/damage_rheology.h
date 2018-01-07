/*
  Copyright (C) 2014 - 2018 by the authors of the ASPECT code.

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


#ifndef __aspect__model_damage_rheology_h
#define __aspect__model_damage_rheology_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/std_cxx1x/array.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional output fields for the dislocation viscosity parameters
     * to be added to the MaterialModel::MaterialModelOutputs structure
     * and filled in the MaterialModel::DamageRheology::evaluate() function.
     */
    template <int dim>
    class DislocationViscosityOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        DislocationViscosityOutputs(const unsigned int n_points);

        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * Dislocation viscosities at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> dislocation_viscosities;

        /**
         * This contains the fraction of the deformation work that is
         * converted to surface energy of grains instead of thermal energy.
         * It is used to reduce the shear heating by this fraction.
         */
        std::vector<double> boundary_area_change_work_fractions;
    };

    namespace Lookup
    {
      /**
       * A base class that can be used to look up material data from an external
       * data source (e.g. a table in a file). The class consists of data members
       * and functions to access this data, but it does not contain the functions
       * to read this data, which has to be implemented in a derived class.
       */
      class MaterialLookup
      {
        public:

          double
          specific_heat(double temperature,
                        double pressure) const;

          double
          density(double temperature,
                  double pressure) const;

          double
          thermal_expansivity(const double temperature,
                              const double pressure) const;

          double
          seismic_Vp(const double temperature,
                     const double pressure) const;

          double
          seismic_Vs(const double temperature,
                     const double pressure) const;

          double
          enthalpy(const double temperature,
                   const double pressure) const;

          double
          dHdT (const double temperature,
                const double pressure) const;

          double
          dHdp (const double temperature,
                const double pressure) const;

          double
          dRhodp (const double temperature,
                  const double pressure) const;

          /**
           * Returns the size of the data tables in pressure (first entry)
           * and temperature (second entry) dimensions.
           */
          std_cxx1x::array<double,2>
          get_pT_steps() const;

        protected:
          /**
           * Access that data value of the property that is stored in table
           * @p values at pressure @p pressure and temperature @p temperature.
           * @p interpol controls whether to perform linear interpolation
           * between the closest data points, or simply use the closest point
           * value.
           */
          double
          value (const double temperature,
                 const double pressure,
                 const dealii::Table<2,
                 double> &values,
                 bool interpol) const;

          /**
           * Find the position in a data table given a temperature.
           */
          double get_nT(double temperature) const;

          /**
           * Find the position in a data table given a pressure.
           */
          double get_np(double pressure) const;

          dealii::Table<2,double> density_values;
          dealii::Table<2,double> thermal_expansivity_values;
          dealii::Table<2,double> specific_heat_values;
          dealii::Table<2,double> vp_values;
          dealii::Table<2,double> vs_values;
          dealii::Table<2,double> enthalpy_values;

          double delta_press;
          double min_press;
          double max_press;
          double delta_temp;
          double min_temp;
          double max_temp;
          unsigned int numtemp;
          unsigned int numpress;
          bool interpolation;
      };

      /**
       * An implementation of the above base class that reads in files created
       * by the HeFESTo software.
       */
      class HeFESToReader : public MaterialLookup
      {
        public:
          HeFESToReader(const std::string &material_filename,
                        const std::string &derivatives_filename,
                        const bool interpol,
                        const MPI_Comm &comm);
      };

      /**
       * An implementation of the above base class that reads in files created
       * by the Perplex software.
       */
      class PerplexReader : public MaterialLookup
      {
        public:
          PerplexReader(const std::string &filename,
                        const bool interpol,
                        const MPI_Comm &comm);
      };
    }

    /**
     * A material model that relies on compositional fields that stand for
     * average grain sizes of a mineral phase and source terms for them that
     * determine the grain size evolution in dependence of the strain rate,
     * temperature, phase transitions, and the creep regime.
     * This material model only works if a compositional field
     * named 'grain_size' is present. In the diffusion
     * creep regime, the viscosity depends on this grain size. We use the grain
     * size evolution laws described in Behn et al., 2009. Implications of grain
     * size evolution on the seismic structure of the oceanic upper mantle, Earth
     * Planet. Sci. Letters, 282, 178–189. Other material parameters are either
     * prescribed similar to the 'simple' material model, or read from data files
     * that were generated by the Perplex or Hefesto software. The material model
     * is described in more detail in "Dannberg et al., 2017. The importance of
     * grain size to mantle dynamics and seismological observations. Geochemistry,
     * Geophysics, Geosystems.", which is the canonical reference for this
     * material model.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class DamageRheology : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the contuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        virtual bool is_compressible () const;

        virtual double reference_viscosity () const;

        virtual void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
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


        virtual
        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;
        /**
         * @}
         */

        /**
         * Returns the enthalpy as calculated by HeFESTo.
         */
        virtual double enthalpy (const double      temperature,
                                 const double      pressure,
                                 const std::vector<double> &compositional_fields,
                                 const Point<dim> &position) const;

        /**
         * Returns the enthalpy derivatives for the evaluate function and
         * postprocessors.
         */
        std_cxx1x::array<std::pair<double, unsigned int>,2>
        enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const;

      protected:
        double reference_rho;
        double reference_T;
        double eta;
        double thermal_alpha;
        double reference_specific_heat;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;

        /**
         * The thermal conductivity.
         */
        double k_value;

        /**
         * Parameters controlling the grain size evolution.
         */
        std::vector<double> grain_growth_activation_energy;
        std::vector<double> grain_growth_activation_volume;
        std::vector<double> grain_growth_rate_constant;
        std::vector<double> grain_growth_exponent;
        std::vector<double> reciprocal_required_strain;
        std::vector<double> recrystallized_grain_size;

        /**
         * Parameters controlling the dynamic grain recrystallization
         * (following paleowattmeter).
         */
        bool use_paleowattmeter;
        std::vector<double> grain_boundary_energy;
        std::vector<double> boundary_area_change_work_fraction;
        std::vector<double> geometric_constant;

        /**
         * Parameters controlling the viscosity.
         */
        double dislocation_viscosity_iteration_threshold;
        unsigned int dislocation_viscosity_iteration_number;
        std::vector<double> dislocation_creep_exponent;
        std::vector<double> dislocation_activation_energy;
        std::vector<double> dislocation_activation_volume;
        std::vector<double> dislocation_creep_prefactor;
        std::vector<double> diffusion_creep_exponent;
        std::vector<double> diffusion_activation_energy;
        std::vector<double> diffusion_activation_volume;
        std::vector<double> diffusion_creep_prefactor;
        std::vector<double> diffusion_creep_grain_size_exponent;

        /**
         * Because of the nonlinear nature of this material model many
         * parameters need to be kept within bounds to ensure stability of the
         * solution. These bounds can be adjusted as input parameters.
         */
        double max_temperature_dependence_of_eta;
        double min_eta;
        double max_eta;
        double min_specific_heat;
        double max_specific_heat;
        double min_thermal_expansivity;
        double max_thermal_expansivity;
        unsigned int max_latent_heat_substeps;
        double min_grain_size;
        double pv_grain_size_scaling;

        /**
         * Whether to advect the real grain size, or the logarithm of the
         * grain size. The logarithm reduces jumps.
         */
        bool advect_log_gransize;


        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double diffusion_viscosity (const double      temperature,
                                            const double      pressure,
                                            const std::vector<double>    &compositional_fields,
                                            const SymmetricTensor<2,dim> &,
                                            const Point<dim> &position) const;

        /**
         * This function calculates the dislocation viscosity. For this purpose
         * we need the dislocation component of the strain rate, which we can
         * only compute by knowing the dislocation viscosity. Therefore, we
         * iteratively solve for the dislocation viscosity and update the
         * dislocation strain rate in each iteration using the new value
         * obtained for the dislocation viscosity. The iteration is started
         * with a dislocation viscosity calculated for the whole strain rate
         * unless a guess for the viscosity is provided, which can reduce the
         * number of iterations significantly.
         */
        virtual double dislocation_viscosity (const double      temperature,
                                              const double      pressure,
                                              const std::vector<double>    &compositional_fields,
                                              const SymmetricTensor<2,dim> &strain_rate,
                                              const Point<dim> &position,
                                              const double viscosity_guess = 0) const;

        /**
         * This function calculates the dislocation viscosity for a given
         * dislocation strain rate.
         */
        double dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                                        const double      pressure,
                                                        const std::vector<double> &,
                                                        const SymmetricTensor<2,dim> &dislocation_strain_rate,
                                                        const Point<dim> &position) const;

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

        /**
         * Returns the p-wave velocity as calculated by HeFESTo.
         */
        virtual double seismic_Vp (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const;

        /**
         * Returns the s-wave velocity as calculated by HeFESTo.
         */
        virtual double seismic_Vs (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const;

        /**
         * Rate of grain size growth (Ostwald ripening) or reduction
         * (due to phase transformations) in dependence on temperature
         * pressure, strain rate, mineral phase and creep regime.
         * We use the grain size evolution laws described in Solomatov
         * and Reese, 2008. Grain size variations in the Earth’s mantle
         * and the evolution of primordial chemical heterogeneities,
         * J. Geophys. Res., 113, B07408.
         */
        virtual
        double
        grain_size_change (const double                  temperature,
                           const double                  pressure,
                           const std::vector<double>    &compositional_fields,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const Tensor<1,dim>          &velocity,
                           const Point<dim>             &position,
                           const unsigned int            phase_index,
                           const int                     crossed_transition) const;

        /**
         * Function that defines the phase transition interface
         * (0 above, 1 below the phase transition).This is done
         * individually for each transition and summed up in the end.
         */
        virtual
        double
        phase_function (const Point<dim> &position,
                        const double temperature,
                        const double pressure,
                        const unsigned int phase) const;

        /**
         * Function that returns the phase for a given
         * position, temperature, pressure and compositional
         * field index.
         */
        virtual
        unsigned int
        get_phase_index (const Point<dim> &position,
                         const double temperature,
                         const double pressure) const;

        /**
         * Function that takes an object in the same format
         * as in.composition as argument and converts the
         * vector that corresponds to the grain size to its
         * logarithms and back and limits the grain size to
         * a global minimum.
         * @in normal_to_log: if true, convert from the grain
         * size to its logarithm, otherwise from log to grain
         * size
         */
        virtual
        void
        convert_log_grain_size (const bool normal_to_log,
                                std::vector<double> &compositional_fields) const;

        /**
         * list of depth, width and Clapeyron slopes for the different phase
         * transitions and in which phase they occur
         */
        std::vector<double> transition_depths;
        std::vector<double> transition_temperatures;
        std::vector<double> transition_slopes;
        std::vector<double> transition_widths;


        /**
         * The following variables are properties of the material files
         * we read in.
         */
        std::string datadirectory;
        std::vector<std::string> material_file_names;
        std::vector<std::string> derivatives_file_names;
        unsigned int n_material_data;
        bool use_table_properties;
        bool use_enthalpy;
        bool use_bilinear_interpolation;


        /**
         * The format of the provided material files. Currently we support
         * the PERPLEX and HeFESTo data formats.
         */
        enum formats
        {
          perplex,
          hefesto
        } material_file_format;

        /**
         * List of pointers to objects that read and process data we get from
         * material data files. There is one pointer/object per compositional
         * field provided.
         */
        std::vector<std_cxx1x::shared_ptr<MaterialModel::Lookup::MaterialLookup> > material_lookup;
    };

  }
}

#endif
