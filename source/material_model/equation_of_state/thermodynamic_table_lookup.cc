/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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
        // This equation of state adds the fraction volume fractions of
        // minerals from different lookups which have the same phase name.
        std::set<std::string> set_phase_volume_column_names;
        n_material_lookups = material_file_names.size();

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
            for (unsigned int j = 0; j < phase_volume_column_names.size(); j++)
              {
                // iterate over the present unique_phase_names object
                // to find phase_volume_column_names[j].
                std::vector<std::string>::iterator it = std::find(unique_phase_names.begin(),
                                                                  unique_phase_names.end(),
                                                                  phase_volume_column_names[j]);

                // If phase_volume_column_names[j] already exists in unique_phase_names,
                // std::distance finds its index. Otherwise, std::distance will return
                // the size of the present object, which is the index where we are
                // about to push the new name. Either way, this is the index we want
                // to add to the unique_phase_indices[i] vector.
                unsigned int i_unique = std::distance(unique_phase_names.begin(), it);
                unique_phase_indices[i].push_back(i_unique);

                // If phase_volume_column_names[j] did not already exist
                // in unique_phase_names, we add it here.
                if (it == unique_phase_names.end())
                  unique_phase_names.push_back(phase_volume_column_names[j]);
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
      double
      ThermodynamicTableLookup<dim>::
      enthalpy (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &/*position*/) const
      {
        double enthalpy = 0.0;

        if (n_material_lookups == 1)
          {
            enthalpy = material_lookup[0]->enthalpy(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_enthalpy = material_lookup[0]->enthalpy(temperature,pressure);
            enthalpy = background_enthalpy;
            for (unsigned int i = 0; i < compositional_fields.size(); ++i)
              enthalpy += compositional_fields[i] *
                          (material_lookup[i+1]->enthalpy(temperature,pressure) - background_enthalpy);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; i++)
              enthalpy += compositional_fields[i] * material_lookup[i]->enthalpy(temperature,pressure);
          }
        return enthalpy;
      }


      template <int dim>
      double
      ThermodynamicTableLookup<dim>::
      seismic_Vp (const double      temperature,
                  const double      pressure,
                  const std::vector<double> &compositional_fields,
                  const Point<dim> &/*position*/) const
      {
        double vp = 0.0;

        if (n_material_lookups == 1)
          {
            vp = material_lookup[0]->seismic_Vp(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_vp = material_lookup[0]->seismic_Vp(temperature,pressure);
            vp = background_vp;
            for (unsigned int i = 0; i < compositional_fields.size(); ++i)
              vp += compositional_fields[i] *
                    (material_lookup[i+1]->seismic_Vp(temperature,pressure) - background_vp);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; i++)
              vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature,pressure);
          }
        return vp;
      }



      template <int dim>
      double
      ThermodynamicTableLookup<dim>::
      seismic_Vs (const double      temperature,
                  const double      pressure,
                  const std::vector<double> &compositional_fields,
                  const Point<dim> &/*position*/) const
      {
        double vs = 0.0;

        if (n_material_lookups == 1)
          {
            vs = material_lookup[0]->seismic_Vs(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_vs = material_lookup[0]->seismic_Vs(temperature,pressure);
            vs = background_vs;
            for (unsigned int i = 0; i < compositional_fields.size(); ++i)
              vs += compositional_fields[i] *
                    (material_lookup[i+1]->seismic_Vs(temperature,pressure) - background_vs);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; i++)
              vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature,pressure);
          }
        return vs;
      }



      template <int dim>
      double
      ThermodynamicTableLookup<dim>::
      density (const double temperature,
               const double pressure,
               const std::vector<double> &compositional_fields, /*composition*/
               const Point<dim> &) const
      {
        double rho = 0.0;
        if (n_material_lookups == 1)
          {
            rho = material_lookup[0]->density(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_density = material_lookup[0]->density(temperature,pressure);
            rho = background_density;
            for (unsigned int i = 0; i < compositional_fields.size(); ++i)
              rho += compositional_fields[i] *
                     (material_lookup[i+1]->density(temperature,pressure) - background_density);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; ++i)
              rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
          }

        return rho;
      }



      template <int dim>
      bool
      ThermodynamicTableLookup<dim>::
      is_compressible () const
      {
        return true;
      }



      template <int dim>
      double
      ThermodynamicTableLookup<dim>::
      compressibility (const double temperature,
                       const double pressure,
                       const std::vector<double> &compositional_fields,
                       const Point<dim> &position) const
      {
        double dRhodp = 0.0;
        if (n_material_lookups == 1)
          {
            dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
            dRhodp = background_dRhodp;
            for (unsigned int i = 0; i < compositional_fields.size(); ++i)
              dRhodp += compositional_fields[i] *
                        (material_lookup[i+1]->dRhodp(temperature,pressure) - background_dRhodp);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; i++)
              dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
          }

        const double rho = density(temperature,pressure,compositional_fields,position);
        return (1/rho)*dRhodp;
      }



      template <int dim>
      double
      ThermodynamicTableLookup<dim>::
      thermal_expansion_coefficient (const double      temperature,
                                     const double      pressure,
                                     const std::vector<double> &compositional_fields,
                                     const Point<dim> &/*position*/) const
      {
        double alpha = 0.0;

        if (n_material_lookups == 1)
          {
            alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
            alpha = background_alpha;
            for (unsigned int i = 0; i<compositional_fields.size(); ++i)
              alpha += compositional_fields[i] *
                       (material_lookup[i+1]->thermal_expansivity(temperature,pressure) - background_alpha);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; ++i)
              alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
          }

        alpha = std::max(std::min(alpha,max_thermal_expansivity),min_thermal_expansivity);
        return alpha;
      }



      template <int dim>
      double
      ThermodynamicTableLookup<dim>::
      specific_heat (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &/*position*/) const
      {
        double cp = 0.0;

        if (n_material_lookups == 1)
          {
            cp = material_lookup[0]->specific_heat(temperature,pressure);
          }
        else if (n_material_lookups == compositional_fields.size() + 1)
          {
            const double background_cp = material_lookup[0]->specific_heat(temperature,pressure);
            cp = background_cp;
            for (unsigned int i = 0; i < compositional_fields.size(); ++i)
              cp += compositional_fields[i] *
                    (material_lookup[i+1]->specific_heat(temperature,pressure) - background_cp);
          }
        else
          {
            for (unsigned i = 0; i < n_material_lookups; ++i)
              cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
          }

        cp = std::max(std::min(cp,max_specific_heat),min_specific_heat);
        return cp;
      }



      template <int dim>
      std::array<std::pair<double, unsigned int>,2>
      ThermodynamicTableLookup<dim>::
      enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const
      {
        std::array<std::pair<double, unsigned int>,2> derivative;

        if (in.current_cell.state() == IteratorState::valid)
          {
            // get the pressures and temperatures at the vertices of the cell
            const QTrapez<dim> quadrature_formula;
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
          }

        return derivative;
      }

      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      fill_phase_volume_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                   NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out) const
      {
        // In the following function,
        // the index i corresponds to the ith compositional field
        // the index j corresponds to the jth phase in the lookup
        // the index k corresponds to the kth evaluation point
        std::vector<std::vector<double>> volume_fractions(unique_phase_names.size(), std::vector<double>(in.n_evaluation_points(), 0.));

        if (n_material_lookups == 1)
          {
            // if there is only one lookup, unique_phase_names is in the same order as the lookup
            for (unsigned int j = 0; j < unique_phase_indices[0].size(); j++)
              {
                for (unsigned int k = 0; k < in.n_evaluation_points(); k++)
                  volume_fractions[j][k] = material_lookup[0]->phase_volume_fraction(j,in.temperature[k],in.pressure[k]);
              }
          }
        else if (n_material_lookups == in.composition[0].size() + 1) // background field
          {
            // The first material lookup corresponds to a background composition
            // We first look up the volume fractions corresponding to this background field
            // (assuming the domain is filled 100% with this single composition)
            std::vector<std::vector<double>> background_volume_fractions(unique_phase_names.size(),
                                                                         std::vector<double>(in.n_evaluation_points(), 0.));

            // We can now loop through the other material models
            for (unsigned int i = 0; i < n_material_lookups; i++)
              {
                for (unsigned int j = 0; j < unique_phase_indices[i].size(); j++)
                  {
                    for (unsigned int k = 0; k < in.n_evaluation_points(); k++)
                      if (i == 0)
                        {
                          background_volume_fractions[unique_phase_indices[0][j]][k] += material_lookup[0]->phase_volume_fraction(j,in.temperature[k],in.pressure[k]);
                          volume_fractions[unique_phase_indices[0][j]][k] = background_volume_fractions[unique_phase_indices[0][j]][k];
                        }
                      else
                        volume_fractions[unique_phase_indices[i][j]][k] += in.composition[k][i-1] * (material_lookup[i]->phase_volume_fraction(j,in.temperature[k],in.pressure[k]) - background_volume_fractions[unique_phase_indices[i][j]][k]);
                  }
              }
          }
        else if (n_material_lookups == in.composition[0].size())
          {
            for (unsigned i = 0; i < n_material_lookups; i++)
              {
                for (unsigned int j = 0; j < unique_phase_indices[i].size(); j++)
                  {
                    for (unsigned int k = 0; k < in.n_evaluation_points(); k++)
                      volume_fractions[unique_phase_indices[i][j]][k] = in.composition[k][i] * material_lookup[i]->phase_volume_fraction(j,in.temperature[k],in.pressure[k]);
                  }
              }
          }
        else
          {
            AssertThrow (false,
                         ExcMessage("The number of material lookups must be equal to "
                                    "one, the number of compositional fields, or the number "
                                    "of compositional fields plus one (if using a background field)."));
          }
        phase_volume_fractions_out->output_values = volume_fractions;
      }

      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      evaluate_enthalpy_dependent_properties(const MaterialModel::MaterialModelInputs<dim> &in,
                                             const unsigned int i,
                                             const double pressure,
                                             const double average_density,
                                             const double average_temperature,
                                             const std::array<std::pair<double, unsigned int>,2> dH,
                                             MaterialModel::MaterialModelOutputs<dim> &out) const
      {

        if (this->get_adiabatic_conditions().is_initialized()
            && (in.current_cell.state() == IteratorState::valid)
            && (dH[0].second > 0)
            && (dH[1].second > 0))
          {
            out.thermal_expansion_coefficients[i] = (1 - average_density * dH[1].first) / average_temperature;
            out.specific_heat[i] = dH[0].first;
          }
        else
          {
            if (n_material_lookups == 1)
              {
                out.thermal_expansion_coefficients[i] = (1 - out.densities[i] * material_lookup[0]->dHdp(in.temperature[i],pressure)) / in.temperature[i];
                out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],pressure);
              }
            else
              {
                ExcNotImplemented();
              }
          }

        out.thermal_expansion_coefficients[i] = std::max(std::min(out.thermal_expansion_coefficients[i],max_thermal_expansivity),min_thermal_expansivity);
        out.specific_heat[i] = std::max(std::min(out.specific_heat[i],max_specific_heat),min_specific_heat);

      }

      template <int dim>
      void
      ThermodynamicTableLookup<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        double average_temperature(0.0);
        double average_density(0.0);
        std::array<std::pair<double, unsigned int>,2> dH;

        if (latent_heat)
          {
            for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
              {
                average_temperature += in.temperature[i];
                average_density += out.densities[i];
              }
            average_temperature /= in.n_evaluation_points();
            average_density /= in.n_evaluation_points();
            dH = enthalpy_derivative(in);
          }

        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            // Use the adiabatic pressure instead of the real one, because of oscillations
            const double pressure = (this->get_adiabatic_conditions().is_initialized())
                                    ?
                                    this->get_adiabatic_conditions().pressure(in.position[i])
                                    :
                                    in.pressure[i];

            out.densities[i] = density(in.temperature[i], pressure, in.composition[i], in.position[i]);
            out.compressibilities[i] = compressibility(in.temperature[i], pressure, in.composition[i], in.position[i]);

            if (latent_heat)
              {
                /* We separate the calculation of specific heat and thermal expansivity,
                 * because they depend on cell-wise averaged values
                 */
                evaluate_enthalpy_dependent_properties(in, i, pressure, average_density, average_temperature, dH, out);
              }
            else
              {
                out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient(in.temperature[i], pressure, in.composition[i], in.position[i]);
                out.specific_heat[i] = specific_heat(in.temperature[i], pressure, in.composition[i], in.position[i]);
              }

            // fill seismic velocities outputs if they exist
            if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim> >())
              {
                seismic_out->vp[i] = seismic_Vp(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
                seismic_out->vs[i] = seismic_Vs(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
              }

            // fill phase volume outputs if they exist
            if (NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out = out.template get_additional_output<NamedAdditionalMaterialOutputs<dim> >())
              fill_phase_volume_fractions(in, phase_volume_fractions_out);
          }
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
                           "be in order with the ordering of the fields). ");


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

        prm.declare_entry ("Minimum specific heat", "500.",
                           Patterns::Double (0.),
                           "The minimum specific heat that is allowed in the whole model domain. "
                           "Units: J/kg/K.");
        prm.declare_entry ("Maximum specific heat", "6000.",
                           Patterns::Double (0.),
                           "The maximum specific heat that is allowed in the whole model domain. "
                           "Units: J/kg/K.");
        prm.declare_entry ("Minimum thermal expansivity", "1e-5",
                           Patterns::Double (),
                           "The minimum thermal expansivity that is allowed in the whole model domain. "
                           "Units: 1/K.");
        prm.declare_entry ("Maximum thermal expansivity", "1e-3",
                           Patterns::Double (),
                           "The maximum thermal expansivity that is allowed in the whole model domain. "
                           "Units: 1/K.");

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
        derivatives_file_names = Utilities::split_string_list(prm.get ("Derivatives file names"));
        use_bilinear_interpolation   = prm.get_bool ("Bilinear interpolation");
        latent_heat                  = prm.get_bool ("Latent heat");

        min_specific_heat            = prm.get_double ("Minimum specific heat");
        max_specific_heat            = prm.get_double ("Maximum specific heat");
        min_thermal_expansivity      = prm.get_double ("Minimum thermal expansivity");
        max_thermal_expansivity      = prm.get_double ("Maximum thermal expansivity");
        max_latent_heat_substeps      = prm.get_integer ("Maximum latent heat substeps");

        if (prm.get ("Material file format") == "perplex")
          material_file_format = perplex;
        else if (prm.get ("Material file format") == "hefesto")
          material_file_format = hefesto;
        else
          AssertThrow (false, ExcNotImplemented());

        // Do some error checking
        AssertThrow ((material_file_names.size() == 1) ||
                     (material_file_names.size() == this->n_compositional_fields()) ||
                     (material_file_names.size() == this->n_compositional_fields() + 1),
                     ExcMessage("This material model expects either one material data file, or as many files as compositional fields, "
                                "or as many files as compositional fields plus one (in which case the first file "
                                "is assumed to contain a background composition). This condition is not fulfilled. You "
                                "prescribed " + Utilities::int_to_string(material_file_names.size()) + " material data files, but there are " +
                                Utilities::int_to_string(this->n_compositional_fields()) + " compositional fields."));

        if (latent_heat)
          AssertThrow (material_file_names.size() == 1,
                       ExcMessage("This formalism is currently only implemented for one material "
                                  "table."));
      }


      template <int dim>
      void
      ThermodynamicTableLookup<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
      {
        if (out.template get_additional_output<NamedAdditionalMaterialOutputs<dim> >() == nullptr)
          {
            const unsigned int n_points = out.n_evaluation_points();
            out.additional_outputs.push_back(
              std_cxx14::make_unique<MaterialModel::NamedAdditionalMaterialOutputs<dim>> (unique_phase_names, n_points));
          }

        if (out.template get_additional_output<SeismicAdditionalOutputs<dim> >() == nullptr)
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
