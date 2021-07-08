/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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
       * An equation of state class that reads thermodynamic properties
       * from pressure-temperature tables in input files. These input files
       * can be created using codes such as Perple_X or HeFESTo.
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
           * Return whether the model is compressible or not.  Incompressibility
           * does not necessarily imply that the density is constant; rather, it
           * may still depend on temperature or pressure. In the current
           * context, compressibility means whether we should solve the continuity
           * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
           * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
           */
          bool is_compressible () const;

          void fill_mass_and_volume_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                               std::vector<std::vector<double>> &mass_fractions,
                                               std::vector<std::vector<double>> &volume_fractions) const;

          /**
          * Function to compute the thermodynamic properties in @p out given the
          * inputs in @p in over all evaluation points.
          * This function also fills the mass_fraction and volume_fraction vectors.
          */
          void
          evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                   std::vector<MaterialModel::EquationOfStateOutputs<dim>> &eos_outputs) const;

          /**
          * Function to fill the seismic velocity and phase volume additional outputs
          */
          void
          fill_additional_outputs(const MaterialModel::MaterialModelInputs<dim> &in,
                                  const std::vector<std::vector<double>> &volume_fractions,
                                  MaterialModel::MaterialModelOutputs<dim> &out) const;

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

          void
          create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const;


        private:
          bool has_background;

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
           * The maximum number of substeps over the temperature pressure range
           * to calculate the averaged enthalpy gradient over a cell
           */
          unsigned int max_latent_heat_substeps;

          /**
           * The format of the provided material files. Currently we support
           * the Perple_X and HeFESTo data formats.
           */
          enum formats
          {
            perplex,
            hefesto
          } material_file_format;

          /**
           * List of pointers to objects that read and process data we get from
           * material data files. There is one pointer/object per lookup file.
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

          /**
           * Compute the specific heat and thermal expansivity using the pressure
           * and temperature derivatives of the specific enthalpy.
           * This evaluation incorporates the effects of latent heat production.
           */
          void evaluate_thermal_enthalpy_derivatives(const MaterialModel::MaterialModelInputs<dim> &in,
                                                     std::vector<MaterialModel::EquationOfStateOutputs<dim>> &eos_outputs) const;

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
          enthalpy_derivatives (const typename Interface<dim>::MaterialModelInputs &in) const;

          void fill_seismic_velocities (const MaterialModel::MaterialModelInputs<dim> &in,
                                        const std::vector<double> &composite_densities,
                                        const std::vector<std::vector<double>> &volume_fractions,
                                        SeismicAdditionalOutputs<dim> *seismic_out) const;

          /**
          * This function uses the MaterialModelInputs &in to fill the output_values
          * of the phase_volume_fractions_out output object with the volume
          * fractions of each of the unique phases at each of the evaluation points.
          * These volume fractions are obtained from the Perple_X-derived
          * pressure-temperature lookup tables.
          * The filled output_values object is a vector of vector<double>;
          * the outer vector is expected to have a size that equals the number
          * of unique phases, the inner vector is expected to have a size that
          * equals the number of evaluation points.
          */
          void fill_phase_volume_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                            const std::vector<std::vector<double>> &volume_fractions,
                                            NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out) const;
      };
    }
  }
}

#endif
