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

#ifndef _aspect_material_model_utilities_h
#define _aspect_material_model_utilities_h

#include <aspect/global.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  template <int dim> class SimulatorAccess;

  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A namespace in which we define utility functions that
     * might be used in many different places in the material
     * model to prevent code duplication.
     */
    namespace MaterialUtilities
    {
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
            specific_heat(const double temperature,
                          const double pressure) const;

            double
            density(const double temperature,
                    const double pressure) const;

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

            /**
             * Computes the derivative of enthalpy for temperature, using the
             * resolution of the read-in table to compute a finite-difference
             * approximation of the derivative.
             */
            double
            dHdT (const double temperature,
                  const double pressure) const;

            /**
            * Computes the derivative of enthalpy for pressure, using the
            * resolution of the read-in table to compute a finite-difference
            * approximation of the derivative.
            */
            double
            dHdp (const double temperature,
                  const double pressure) const;

            /**
             * Compute the enthalpy derivatives for temperature and pressure
             * given a set of temperature and pressure points, which will be
             * used as support points for the finite difference scheme. This
             * is useful to not 'miss' phase transitions that are not resolved in
             * the dHdT and dHdp functions. The third argument represents
             * the number of substeps taken to compute this average. A number
             * larger than one means the temperature-pressure range that is spanned
             * by the first two input arguments is seperated into @p n_substeps
             * equally spaced pressure-temperature steps, the derivatives are
             * computed for each substep and then averaged.
             */
            std::array<std::pair<double, unsigned int>,2>
            enthalpy_derivatives(const std::vector<double> &temperatures,
                                 const std::vector<double> &pressures,
                                 const unsigned int n_substeps = 1) const;

            double
            dRhodp (const double temperature,
                    const double pressure) const;

            /**
             * Returns the size of the data tables in pressure (first entry)
             * and temperature (second entry) dimensions.
             */
            std::array<double,2>
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
                   const Table<2, double> &values,
                   const bool interpol) const;

            /**
             * Find the position in a data table given a temperature.
             */
            double get_nT(const double temperature) const;

            /**
             * Find the position in a data table given a pressure.
             */
            double get_np(const double pressure) const;

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
            unsigned int n_temperature;
            unsigned int n_pressure;
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
       * For multicomponent material models: Given a vector of of compositional
       * fields of length N, this function returns a vector of volume fractions
       * of length N+1, corresponding to the volume fraction of a ``background
       * material'' as the first entry, and volume fractions for each of the input
       * fields as the following entries. The returned vector will sum to one.
       * If the sum of the compositional_fields is greater than
       * one, we assume that there is no background mantle (i.e., that field value
       * is zero). Otherwise, the difference between the sum of the compositional
       * fields and 1.0 is assumed to be the amount of background mantle.
       * Optionally, one can input a component mask that determines which of the
       * compositional fields to use during the computation (e.g. because
       * some fields contain non-volumetric quantities like strain,
       * porosity, or trace elements). By default, all fields are included.
       */
      std::vector<double>
      compute_volume_fractions(const std::vector<double> &compositional_fields,
                               const ComponentMask &field_mask = ComponentMask());



      /**
       * For multicomponent material models:
       * Enumeration for selecting which averaging scheme to use when
       * averaging the properties of different compositional fields.
       * Select between harmonic, arithmetic, geometric, and
       * maximum_composition. The max composition scheme simply uses the
       * viscosity of whichever field has the highest volume fraction.
       */
      enum CompositionalAveragingOperation
      {
        harmonic,
        arithmetic,
        geometric,
        maximum_composition
      };



      /**
       * Read the compositional averaging operation from the parameter file,
       * using the parameter name given in @p parameter_name, and return the
       * enum that corresponds to this operation.
       */
      CompositionalAveragingOperation
      parse_compositional_averaging_operation (const std::string &parameter_name,
                                               const ParameterHandler &prm);



      /**
       * For multicomponent material models:
       * Material models compute output quantities such as the viscosity, the
       * density, etc. For some models, these values depend strongly on the
       * composition, and more than one compositional field might have nonzero
       * values at a given quadrature point. This means that properties have to
       * be averaged based on the fractions of each compositional field present.
       * This function performs this type of averaging. The averaging is based
       * on the choice in @p average_type. Averaging is conducted over the
       * compositional fields given in @p volume_fractions. This means that
       * @p volume_fractions and @p parameter_values have to have the same size,
       * which would typically be the number of compositional fields used in the
       * simulation (with the potential addition of a background field, in case
       * the composition does not add up to 1). However, one might not want to
       * average over all fields, as in some cases compositional fields do not
       * represent a rock type, but other tracked quantities like the finite
       * strain, so the implementation is independent of the number of entries in
       * @p volume_fractions.
       */
      double average_value (const std::vector<double> &volume_fractions,
                            const std::vector<double> &parameter_values,
                            const CompositionalAveragingOperation &average_type);



      /**
       * A data structure with all inputs for the
       * MaterialModel::Interface::compute_drucker_prager_yielding() method.
       */
      struct DruckerPragerInputs
      {
        /**
         * Constructor. Initializes the various variables of this structure
         * with the input values. By default, there is no maximum yield
         * strength, so the parameter is set to infinity.
         */
        DruckerPragerInputs(const double cohesion,
                            const double friction_angle,
                            const double pressure,
                            const double effective_strain_rate,
                            const double max_yield_strength = std::numeric_limits<double>::infinity());

        const double cohesion;
        const double friction_angle;
        const double pressure;
        const double effective_strain_rate;
        const double max_yield_strength;
      };



      /**
       * A data structure with all outputs computed by the
       * MaterialModel::Interface::compute_drucker_prager_yielding() method.
       */
      struct DruckerPragerOutputs
      {
        /**
         * Constructor. Initializes the various variables of this structure to
         * NaNs.
         */
        DruckerPragerOutputs();

        double yield_strength;
        double plastic_viscosity;
        double viscosity_pressure_derivative;
      };



      /**
       * For material models with plasticity:
       * Function to compute the material properties in @p out given the
       * inputs in @p in according to the the Drucker-Prager yield criterion.
       */
      template <int dim>
      void
      compute_drucker_prager_yielding (const DruckerPragerInputs &in,
                                       DruckerPragerOutputs &out);



      /**
       * A data structure with all inputs for the
       * MaterialModel::Interface::compute_phase_function() method.
       */
      template <int dim>
      struct PhaseFunctionInputs
      {
        /**
        * Constructor. Initializes the various variables of this
        * structure with the input values.
        */
        PhaseFunctionInputs(const double temperature,
                            const double pressure,
                            const Point<dim> &position,
                            const unsigned int phase_index,
                            const unsigned int composition_index);

        double temperature;
        double pressure;
        Point<dim> position;
        unsigned int phase_index;
        unsigned int composition_index;
      };

      /**
             * A simplified, incompressible equation of state where the density depends linearly
             * on temperature and composition, using the equation
             * $\rho(p,T,\mathfrak c) = \left(1-\alpha (T-T_0)\right)\rho_0 + \sum_i \Delta\rho_i \; \mathfrak c_i.$ "
             * There is no pressure-dependence of the density, and all other material properties
             * relating to the equation of state are assumed to be constant and identical for each
             * composition.
             */
      template <int dim>
      class PhaseFunction: public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * For material models with phase changes:
           * Function to compute the percentage of material that has
           * already undergone the phase transition to the higher-pressure
           * material. This is done individually for each transition and
           * summed up in the end.
           */
          double phase_function (const PhaseFunctionInputs<dim> &in) const;

          /**
           * For material models with phase changes:
           * Function to compute the derivative of the phase function.
           */
          double phase_function_derivative (const PhaseFunctionInputs<dim> &in) const;

          unsigned int n_phase_transitions () const;

          double get_transition_slope (const unsigned int phase_index) const;

          /**
           * Declare the parameters this class takes through input files.
           * The optional parameter @p n_compositions determines the maximum
           * number of compositions the equation of state is set up with,
           * in other words, how many compositional fields influence the
           * density.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * The optional parameter @p n_compositions determines the maximum
           * number of compositions the equation of state is set up with,
           * and should have the same value as the parameter with the same
           * name in the declare_parameters() function.
           */
          void
          parse_parameters (ParameterHandler &prm);


        private:
          /**
           * A function to compute the phase transition pressure and pressure width
           * if the phase transitions are defined by the user in terms of depth.
           * This function should be moved to the material model utilities.
           */
          virtual
          std::pair<double, double>
          transition_depth_to_pressure (const Point<dim> &position,
                                        const int phase) const;

          // list of depth (or pressure), width and Clapeyron slopes
          // for the different phase transitions
          std::vector<double> transition_depths;
          std::vector<double> transition_pressures;
          std::vector<double> transition_temperatures;
          std::vector<double> transition_widths;
          std::vector<double> transition_pressure_widths;
          std::vector<double> transition_slopes;

          bool use_depth_instead_of_pressure;
      };
    }
  }
}


#endif
