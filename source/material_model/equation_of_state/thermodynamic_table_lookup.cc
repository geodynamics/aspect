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


#include <aspect/material_model/equation_of_state/thermodynamic_table_lookup.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      template <int dim>
      void
      ThermodynamicTableLookup<dim>::initialize()
      {
        // This model allows the user to provide several PerpleX or HeFESTo
        // P-T lookup files, each of which corresponds to a different material.
        // The thermodynamic properties and stable mineral phases within each material
        // will in general be different, and must be averaged in a reasonable way.

        // In the following, we populate the material_lookups.
        // We also populate the unique_phase_names vector, which is used to
        // calculate the phase volumes if they are provided by the lookup.
        // This equation of state sums the volume fractions of
        // phases from different lookups if they have the same phase name.
        std::set<std::string> set_phase_volume_column_names;

        // Resize the unique_phase_indices object
        unique_phase_indices.resize(n_material_lookups, std::vector<unsigned int>());

        for (unsigned i = 0; i < n_material_lookups; i++)
          {
            if (material_file_format == perplex)
              material_lookup
              .push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>(data_directory+material_file_names[i],
                         use_bilinear_interpolation,
                         this->get_mpi_communicator()));
            else if (material_file_format == hefesto)
              material_lookup
              .push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::HeFESToReader>(data_directory+material_file_names[i],
                         data_directory+derivatives_file_names[i],
                         use_bilinear_interpolation,
                         this->get_mpi_communicator()));
            else
              AssertThrow (false, ExcNotImplemented());

            // Here we look up all of the column names and insert
            // unique names into the unique_phase_names vector and
            // filling the unique_phase_indices object.
            std::vector<std::string> phase_volume_column_names = material_lookup[i]->phase_volume_column_names();
            for (const auto &phase_volume_column_name : phase_volume_column_names)
              {
                // iterate over the present unique_phase_names object
                // to find phase_volume_column_name.
                std::vector<std::string>::iterator it = std::find(unique_phase_names.begin(),
                                                                  unique_phase_names.end(),
                                                                  phase_volume_column_name);

                // If phase_volume_column_name already exists in unique_phase_names,
                // std::distance finds its index. Otherwise, std::distance will return
                // the size of the present object, which is the index where we are
                // about to push the new name. Either way, this is the index we want
                // to add to the unique_phase_indices[i] vector.
                unsigned int i_unique = std::distance(unique_phase_names.begin(), it);
                unique_phase_indices[i].push_back(i_unique);

                // If phase_volume_column_name did not already exist
                // in unique_phase_names, we add it here.
                if (it == unique_phase_names.end())
                  unique_phase_names.push_back(phase_volume_column_name);
              }
          }
      }



      template <int dim>
      unsigned int
      ThermodynamicTableLookup<dim>::number_of_lookups() const
      {
        return n_material_lookups;
      }



      template <int dim>
      bool
      ThermodynamicTableLookup<dim>::
      is_compressible () const
      {
        return true;
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      fill_seismic_velocities (const MaterialModel::MaterialModelInputs<dim> &in,
                               const std::vector<double> &composite_densities,
                               const std::vector<std::vector<double>> &volume_fractions,
                               SeismicAdditionalOutputs<dim> *seismic_out) const
      {
        // This function returns the Voigt-Reuss-Hill averages of the
        // seismic velocities of the different materials.

        // Now we calculate the averaged moduli.
        // mu = rho*Vs^2; K_s = rho*Vp^2 - 4./3.*mu
        // The Voigt average is an arithmetic volumetric average,
        // while the Reuss average is a harmonic volumetric average.
        for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
          {
            if (material_lookup.size() == 1)
              {
                seismic_out->vs[i] = material_lookup[0]->seismic_Vs(in.temperature[i],in.pressure[i]);
                seismic_out->vp[i] = material_lookup[0]->seismic_Vp(in.temperature[i],in.pressure[i]);
              }
            else
              {
                double k_voigt = 0.;
                double mu_voigt = 0.;
                double invk_reuss = 0.;
                double invmu_reuss = 0.;

                for (unsigned int j = 0; j < material_lookup.size(); ++j)
                  {
                    const double mu = material_lookup[j]->density(in.temperature[i],in.pressure[i])*std::pow(material_lookup[j]->seismic_Vs(in.temperature[i],in.pressure[i]), 2.);
                    const double k =  material_lookup[j]->density(in.temperature[i],in.pressure[i])*std::pow(material_lookup[j]->seismic_Vp(in.temperature[i],in.pressure[i]), 2.) - 4./3.*mu;

                    k_voigt += volume_fractions[i][j] * k;
                    mu_voigt += volume_fractions[i][j] * mu;
                    invk_reuss += volume_fractions[i][j] / k;
                    invmu_reuss += volume_fractions[i][j] / mu;
                  }

                const double k_VRH = (k_voigt + 1./invk_reuss)/2.;
                const double mu_VRH = (mu_voigt + 1./invmu_reuss)/2.;
                seismic_out->vp[i] = std::sqrt((k_VRH + 4./3.*mu_VRH)/composite_densities[i]);
                seismic_out->vs[i] = std::sqrt(mu_VRH/composite_densities[i]);
              }
          }
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      fill_phase_volume_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                   const std::vector<std::vector<double>> &volume_fractions,
                                   NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out) const
      {
        // Each call to material_lookup[j]->phase_volume_fraction(k,temperature,pressure)
        // returns the volume fraction of the kth phase which is present in that material lookup
        // at the requested temperature and pressure.
        // The total volume fraction of each phase at each evaluation point is equal to
        // sum_j (volume_fraction_of_material_j * phase_volume_fraction_in_material_j).
        // In the following function,
        // the index i corresponds to the ith evaluation point
        // the index j corresponds to the jth compositional field
        // the index k corresponds to the kth phase in the lookup
        std::vector<std::vector<double>> phase_volume_fractions(unique_phase_names.size(),
                                                                std::vector<double>(in.n_evaluation_points(), 0.));
        for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
          for (unsigned j = 0; j < material_lookup.size(); ++j)
            for (unsigned int k = 0; k < unique_phase_indices[j].size(); ++k)
              phase_volume_fractions[unique_phase_indices[j][k]][i] += volume_fractions[i][j] * material_lookup[j]->phase_volume_fraction(k,in.temperature[i],in.pressure[i]);

        phase_volume_fractions_out->output_values = phase_volume_fractions;
      }



      template <int dim>
      std::array<std::pair<double, unsigned int>,2>
      ThermodynamicTableLookup<dim>::
      enthalpy_derivatives (const typename Interface<dim>::MaterialModelInputs &in) const
      {
        std::array<std::pair<double, unsigned int>,2> derivative;

        // get the pressures and temperatures at the vertices of the cell
#if DEAL_II_VERSION_GTE(9,3,0)
        const QTrapezoid<dim> quadrature_formula;
#else
        const QTrapez<dim> quadrature_formula;
#endif

        const unsigned int n_q_points = quadrature_formula.size();
        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 quadrature_formula,
                                 update_values);

        std::vector<double> temperatures(n_q_points), pressures(n_q_points);
        fe_values.reinit (in.current_cell);

        fe_values[this->introspection().extractors.temperature]
        .get_function_values (this->get_current_linearization_point(), temperatures);
        fe_values[this->introspection().extractors.pressure]
        .get_function_values (this->get_current_linearization_point(), pressures);

        AssertThrow (n_material_lookups == 1,
                     ExcMessage("This formalism is only implemented for one material "
                                "table."));

        // We have to take into account here that the p,T spacing of the table of material properties
        // we use might be on a finer grid than our model. Because of that we compute the enthalpy
        // derivatives by using finite differences that average over the whole temperature and
        // pressure range that is used in this cell. This way we should not miss any phase transformation.
        derivative = material_lookup[0]->enthalpy_derivatives(temperatures,
                                                              pressures,
                                                              max_latent_heat_substeps);

        return derivative;
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               std::vector<MaterialModel::EquationOfStateOutputs<dim>> &eos_outputs) const
      {
        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          {
            const double pressure = in.pressure[i];
            const double temperature = in.temperature[i];

            for (unsigned int j=0; j<eos_outputs[i].densities.size(); ++j)
              {
                eos_outputs[i].densities[j] = material_lookup[j]->density(temperature, pressure);
                eos_outputs[i].compressibilities[j] = material_lookup[j]->dRhodp(temperature, pressure)/eos_outputs[i].densities[j];

                // Only calculate the non-reactive specific heat and
                // thermal expansivity if latent heat is to be ignored.
                if (!latent_heat)
                  {
                    eos_outputs[i].thermal_expansion_coefficients[j] = material_lookup[j]->thermal_expansivity(temperature, pressure);
                    eos_outputs[i].specific_heat_capacities[j] = material_lookup[j]->specific_heat(temperature, pressure);
                  }

                eos_outputs[i].entropy_derivative_pressure[j] = 0.;
                eos_outputs[i].entropy_derivative_temperature[j] = 0.;
              }
          }
        // Calculate the specific heat and thermal expansivity
        // including the effects of reaction.
        if (latent_heat)
          {
            evaluate_thermal_enthalpy_derivatives(in, eos_outputs);
          }
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      evaluate_thermal_enthalpy_derivatives(const MaterialModel::MaterialModelInputs<dim> &in,
                                            std::vector<MaterialModel::EquationOfStateOutputs<dim>> &eos_outputs) const
      {
        // The second derivatives of the thermodynamic potentials (compressibility, thermal expansivity, specific heat)
        // are dependent not only on the phases present in the assemblage at the given temperature and pressure,
        // but also on any reactions between phases in the assemblage. PerpleX and HeFESTo output only "static" properties
        // (properties not including any reaction effects), and so do not capture the latent heat of reaction.

        // In this material model, we always use a compressibility which includes the effects of reaction,
        // but we allow the user the option to switch on or off thermal (latent heat) effects.
        // If the latent_heat bool is set to true, thermal expansivity and specific heat are calculated from
        // the change in enthalpy with pressure and temperature.

        // There are alternative ways to capture the latent heat effect (by preprocessing the P-T table, for example),
        // which may be a more appropriate approach in some cases, but the latent heat should always be considered if
        // thermodynamic self-consistency is intended.

        // The affected properties are computed using the partial derivatives of the enthalpy
        // with respect to pressure and temperature:
        // thermal expansivity = (1 - rho * (dH/dp)_T) / T
        // specific heat capacity = (dH/dT)_P
        // where the subscript indicates the natural variable which is held constant.
        double average_temperature(0.0);
        double average_density(0.0);
        std::array<std::pair<double, unsigned int>,2> dH;

        for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
          {
            average_temperature += in.temperature[i];
            average_density += eos_outputs[i].densities[0];
          }
        average_temperature /= in.n_evaluation_points();
        average_density /= in.n_evaluation_points();

        if (in.current_cell.state() == IteratorState::valid)
          dH = enthalpy_derivatives(in);

        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          {
            // Use the adiabatic pressure instead of the real one,
            // to stabilize against pressure oscillations in phase transitions
            const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]);

            if ((in.current_cell.state() == IteratorState::valid)
                && (dH[0].second > 0) && (dH[1].second > 0))
              {
                eos_outputs[i].thermal_expansion_coefficients[0] = (1 - average_density * dH[1].first) / average_temperature;
                eos_outputs[i].specific_heat_capacities[0] = dH[0].first;
              }
            else
              {
                eos_outputs[i].thermal_expansion_coefficients[0] = (1 - eos_outputs[i].densities[0] * material_lookup[0]->dHdp(in.temperature[i],pressure)) / in.temperature[i];
                eos_outputs[i].specific_heat_capacities[0] = material_lookup[0]->dHdT(in.temperature[i],pressure);
              }
          }
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::fill_additional_outputs (const MaterialModel::MaterialModelInputs<dim> &in,
                                                              const std::vector<std::vector<double>> &volume_fractions,
                                                              MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        // fill seismic velocity outputs if they exist
        if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim>>())
          fill_seismic_velocities(in, out.densities, volume_fractions, seismic_out);

        // fill phase volume outputs if they exist
        if (NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out = out.template get_additional_output<NamedAdditionalMaterialOutputs<dim>>())
          fill_phase_volume_fractions(in, volume_fractions, phase_volume_fractions_out);
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                           Patterns::DirectoryName (),
                           "The path to the model data. The path may also include the special "
                           "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the `data/' subdirectory of ASPECT. ");
        prm.declare_entry ("Material file names", "pyr-ringwood88.txt",
                           Patterns::List (Patterns::Anything()),
                           "The file names of the material data (material "
                           "data is assumed to be in order with the ordering "
                           "of the compositional fields). Note that there are "
                           "three options on how many files need to be listed "
                           "here: 1. If only one file is provided, it is used "
                           "for the whole model domain, and compositional fields "
                           "are ignored. 2. If there is one more file name than the "
                           "number of compositional fields, then the first file is "
                           "assumed to define a `background composition' that is "
                           "modified by the compositional fields. If there are "
                           "exactly as many files as compositional fields, the fields are "
                           "assumed to represent the fractions of different materials "
                           "and the average property is computed as a sum of "
                           "the value of the compositional field times the "
                           "material property of that field.");
        prm.declare_entry ("Derivatives file names", "",
                           Patterns::List (Patterns::Anything()),
                           "The file names of the enthalpy derivatives data. "
                           "List with as many components as active "
                           "compositional fields (material data is assumed to "
                           "be in order with the ordering of the fields).");
        prm.declare_entry ("Material file format", "perplex",
                           Patterns::Selection ("perplex|hefesto"),
                           "The material file format to be read in the property "
                           "tables.");
        prm.declare_entry ("Bilinear interpolation", "true",
                           Patterns::Bool (),
                           "Whether to use bilinear interpolation to compute "
                           "material properties (slower but more accurate). ");
        prm.declare_entry ("Latent heat", "false",
                           Patterns::Bool (),
                           "Whether to include latent heat effects in the "
                           "calculation of thermal expansivity and specific heat. "
                           "If true, ASPECT follows the approach of Nakagawa et al. 2009, "
                           "using pressure and temperature derivatives of the enthalpy "
                           "to calculate the thermal expansivity and specific heat. "
                           "If false, ASPECT uses the thermal expansivity and "
                           "specific heat values from the material properties table.");
        prm.declare_entry ("Maximum latent heat substeps", "1",
                           Patterns::Integer (1),
                           "The maximum number of substeps over the temperature pressure range "
                           "to calculate the averaged enthalpy gradient over a cell.");
      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::parse_parameters (ParameterHandler &prm)
      {
        data_directory               = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
        material_file_names          = Utilities::split_string_list(prm.get ("Material file names"));
        n_material_lookups           = material_file_names.size();

        derivatives_file_names       = Utilities::split_string_list(prm.get ("Derivatives file names"));
        use_bilinear_interpolation   = prm.get_bool ("Bilinear interpolation");
        latent_heat                  = prm.get_bool ("Latent heat");
        max_latent_heat_substeps     = prm.get_integer ("Maximum latent heat substeps");

        if (prm.get ("Material file format") == "perplex")
          material_file_format       = perplex;
        else if (prm.get ("Material file format") == "hefesto")
          material_file_format       = hefesto;
        else
          AssertThrow (false, ExcNotImplemented());

        if (latent_heat)
          AssertThrow (n_material_lookups == 1,
                       ExcMessage("Isochemical latent heat calculations are only implemented for a single material lookup."));

      }



      template <int dim>
      void
      ThermodynamicTableLookup<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<NamedAdditionalMaterialOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<MaterialModel::NamedAdditionalMaterialOutputs<dim>> (unique_phase_names, n_points));
          }

        if (out.template get_additional_output<SeismicAdditionalOutputs<dim>>() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
#define INSTANTIATE(dim) \
  template class ThermodynamicTableLookup<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
