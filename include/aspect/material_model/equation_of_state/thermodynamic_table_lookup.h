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

#ifndef _aspect_material_model_equation_of_state_thermodynamic_table_lookup_h
#define _aspect_material_model_equation_of_state_thermodynamic_table_lookup_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      using namespace dealii;

      /**
       *
       */
      template <int dim>
      class ThermodynamicTableLookup : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Initialization function. Loads the material data and sets up
           * pointers.
           */
          void
          initialize ();

          /**
           * Returns the number of lookups
           */
          virtual unsigned int number_of_lookups() const;

          /**
           * @name Physical parameters used in the basic equations
           * @{
           */

          virtual double enthalpy (const double temperature,
                                   const double pressure,
                                   const std::vector<double> &compositional_fields,
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

          virtual double seismic_Vp (const double      temperature,
                                     const double      pressure,
                                     const std::vector<double> &compositional_fields,
                                     const Point<dim> &position) const;

          virtual double seismic_Vs (const double      temperature,
                                     const double      pressure,
                                     const std::vector<double> &compositional_fields,
                                     const Point<dim> &position) const;

          void evaluate_enthalpy_dependent_properties(const MaterialModel::MaterialModelInputs<dim> &in,
                                                      const unsigned int i,
                                                      const double pressure,
                                                      const double average_density,
                                                      const double average_temperature,
                                                      const std::array<std::pair<double, unsigned int>,2> dH,
                                                      MaterialModel::MaterialModelOutputs<dim> &out) const;

          /**
          * This function uses the MaterialModelInputs &in to fill the output_values
          * of the phase_volume_fractions_out output object with the volume
          * fractions of each of the unique phases at each of the evaluation points.
          * These volume fractions are obtained from the PerpleX-derived
          * pressure-temperature lookup tables.
          * The filled output_values object is a vector of vector<double>;
          * the outer vector is expected to have a size that equals the number
          * of unique phases, the inner vector is expected to have a size that
          * equals the number of evaluation points.
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
          bool is_compressible () const;
          /**
           * @}
           */

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
                   MaterialModel::MaterialModelOutputs<dim> &out) const;

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
          parse_parameters (ParameterHandler &prm);
          /**
           * @}
           */

          void
          create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


        private:
          unsigned int n_material_lookups;
          bool use_bilinear_interpolation;
          bool latent_heat;

          /**
           * Information about the location of data files.
           */
          std::string data_directory;
          std::vector<std::string> material_file_names;
          std::vector<std::string> derivatives_file_names;

          /**
           * Because of the nonlinear nature of this material model many
           * parameters need to be kept within bounds to ensure stability of the
           * solution. These bounds can be adjusted as input parameters.
           */
          double min_specific_heat;
          double max_specific_heat;
          double min_thermal_expansivity;
          double max_thermal_expansivity;
          unsigned int max_latent_heat_substeps;

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
          std::vector<std::unique_ptr<MaterialModel::MaterialUtilities::Lookup::MaterialLookup> > material_lookup;

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

      };
    }
  }
}

#endif
