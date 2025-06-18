/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/steinberger.h>
#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/thermal_conductivity/constant.h>
#include <aspect/material_model/thermal_conductivity/tosi_stackhouse.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>
#include <aspect/lateral_averaging.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/table.h>
#include <memory>


namespace aspect
{
  namespace MaterialModel
  {
    namespace internal
    {
      LateralViscosityLookup::LateralViscosityLookup(const std::string &filename,
                                                     const MPI_Comm comm)
      {
        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        std::getline(in, temp); // eat first line

        min_depth=1e20;
        max_depth=-1;

        while (!in.eof())
          {
            double visc, depth;
            in >> visc;
            if (in.eof())
              break;
            in >> depth;
            depth *=1000.0;
            std::getline(in, temp);

            min_depth = std::min(depth, min_depth);
            max_depth = std::max(depth, max_depth);

            values.push_back(visc);
          }
        delta_depth = (max_depth-min_depth)/(values.size()-1);
      }

      double LateralViscosityLookup::lateral_viscosity(double depth) const
      {
        depth=std::max(min_depth, depth);
        depth=std::min(depth, max_depth);

        Assert(depth>=min_depth, ExcMessage("ASPECT found a depth less than min_depth."));
        Assert(depth<=max_depth, ExcMessage("ASPECT found a depth greater than max_depth."));
        const unsigned int idx = static_cast<unsigned int>((depth-min_depth)/delta_depth);
        Assert(idx<values.size(), ExcMessage("Attempting to look up a depth with an index that would be out of range. (depth-min_depth)/delta_depth too large."));
        return values[idx];
      }

      int LateralViscosityLookup::get_nslices() const
      {
        return values.size();
      }

      RadialViscosityLookup::RadialViscosityLookup(const std::string &filename,
                                                   const MPI_Comm comm)
      {
        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        min_depth=1e20;
        max_depth=-1;

        while (!in.eof())
          {
            double visc, depth;
            in >> visc;
            if (in.eof())
              break;
            in >> depth;
            depth *=1000.0;
            std::getline(in, temp);

            min_depth = std::min(depth, min_depth);
            max_depth = std::max(depth, max_depth);

            values.push_back(visc);
          }
        delta_depth = (max_depth-min_depth)/(values.size()-1);
      }

      double RadialViscosityLookup::radial_viscosity(double depth) const
      {
        depth=std::max(min_depth, depth);
        depth=std::min(depth, max_depth);

        Assert(depth>=min_depth, ExcMessage("ASPECT found a depth less than min_depth."));
        Assert(depth<=max_depth, ExcMessage("ASPECT found a depth greater than max_depth."));
        const unsigned int idx = static_cast<unsigned int>((depth-min_depth)/delta_depth);
        Assert(idx<values.size(), ExcMessage("Attempting to look up a depth with an index that would be out of range. (depth-min_depth)/delta_depth too large."));
        return values[idx];
      }
    }



    template <int dim>
    void
    Steinberger<dim>::initialize()
    {
      equation_of_state.initialize();

      lateral_viscosity_lookup
        = std::make_unique<internal::LateralViscosityLookup>(data_directory+lateral_viscosity_file_name,
                                                             this->get_mpi_communicator());
      radial_viscosity_lookup
        = std::make_unique<internal::RadialViscosityLookup>(data_directory+radial_viscosity_file_name,
                                                            this->get_mpi_communicator());
      average_temperature.resize(n_lateral_slices);
    }



    template <int dim>
    void
    Steinberger<dim>::
    update()
    {
      if (use_lateral_average_temperature)
        {
          this->get_lateral_averaging().get_temperature_averages(average_temperature);
          for (double temperature : average_temperature)
            AssertThrow(numbers::is_finite(temperature),
                        ExcMessage("In computing depth averages, there is at"
                                   " least one depth band that does not have"
                                   " any quadrature points in it."
                                   " Consider reducing number of depth layers"
                                   " for averaging specified in the parameter"
                                   " file.(Number lateral average bands)"));
        }
    }



    template <int dim>
    double
    Steinberger<dim>::
    viscosity (const double temperature,
               const double /*pressure*/,
               const std::vector<double> &volume_fractions,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &position,
               bool output_plasticity,
               MaterialModel::MaterialModelOutputs<dim> &out,
               const unsigned int i) const
    {
      const double depth = this->get_geometry_model().depth(position);
      const double adiabatic_temperature = this->get_adiabatic_conditions().temperature(position);

      double delta_temperature;
      if (use_lateral_average_temperature)
        {
          const unsigned int idx = static_cast<unsigned int>((average_temperature.size()-1) * depth / this->get_geometry_model().maximal_depth());
          delta_temperature = temperature-average_temperature[idx];
        }
      else
        delta_temperature = temperature-adiabatic_temperature;

      // For an explanation on this formula see the Steinberger & Calderwood 2006 paper
      // We here compute the lateral variation of viscosity due to temperature (thermal_prefactor) as
      // V_lT = exp [-(H/nR)*dT/(T_adiabatic*(T_adiabatic + dT))] as in Eq. 6 of the paper.
      // We get H/nR from the lateral_viscosity_lookup->lateral_viscosity function.
      const double log_thermal_prefactor = -1.0 * lateral_viscosity_lookup->lateral_viscosity(depth) * delta_temperature / (temperature * adiabatic_temperature);

      // Limit the lateral viscosity variation to a reasonable interval
      const double thermal_prefactor = std::max(std::min(std::exp(log_thermal_prefactor), max_lateral_eta_variation), 1/max_lateral_eta_variation);

      const double compositional_prefactor = MaterialUtilities::average_value (volume_fractions, viscosity_prefactors, viscosity_averaging_scheme);

      // Visc_rT = exp[(H/nR)/T_adiabatic], Eq. 7 of the paper
      const double eta_ref = radial_viscosity_lookup->radial_viscosity(depth);

      // Effective viscosity without Drucker-Prager yielding
      // It is computed as radial viscosity profile * thermal prefactors * compositional prefactor
      double effective_viscosity = thermal_prefactor * compositional_prefactor * eta_ref;
      const double pressure = this->get_adiabatic_conditions().pressure(position);

      // Druger-Pager rheology
      if (enable_drucker_prager_rheology)
        {
          // This should be the same as the "strain_rate" output in the visualization postpocessor,
          // which also known as second strain rate invariant or effective deviatoric strain rate.
          const double strain_rate_effective = std::sqrt(std::fabs(second_invariant(deviator(strain_rate))));
          // Calculate non-yielding (viscous) stress magnitude. It should be identical to the
          // "stress second invariant" in the visualization postprocessor.
          const double non_yielding_stress = 2. * effective_viscosity * strain_rate_effective;

          // In the steinberger material model, viscosity does not depend on composition or phases,
          // so we set the compositional index for the Drucker-Prager parameters to 0 and not using
          // phase averaging.
          const Rheology::DruckerPragerParameters drucker_prager_parameters = drucker_prager_plasticity.compute_drucker_prager_parameters(0);

          const double yield_stress = drucker_prager_plasticity.compute_yield_stress(pressure,drucker_prager_parameters);

          // Apply plastic yielding:
          // If the non-yielding stress is greater than the yield stress,
          // rescale the viscosity back to yield surface
          if (non_yielding_stress >= yield_stress)
            {
              effective_viscosity = drucker_prager_plasticity.compute_viscosity(pressure,
                                                                                strain_rate_effective,
                                                                                drucker_prager_parameters,
                                                                                effective_viscosity);
            }

          const std::shared_ptr<PlasticAdditionalOutputs<dim>> plastic_out
            = out.template get_additional_output_object<PlasticAdditionalOutputs<dim>>();

          if (plastic_out != nullptr && output_plasticity)
            {
              plastic_out->cohesions[i] = drucker_prager_parameters.cohesion;
              plastic_out->friction_angles[i] = drucker_prager_parameters.angle_internal_friction;
              plastic_out->yield_stresses[i] = yield_stress;
              plastic_out->yielding[i] = non_yielding_stress >= yield_stress ? 1 : 0;
            }

        }

      return effective_viscosity;
    }



    template <int dim>
    bool
    Steinberger<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible();
    }



    template <int dim>
    void
    Steinberger<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      std::vector<EquationOfStateOutputs<dim>> eos_outputs (in.n_evaluation_points(), equation_of_state.number_of_lookups());
      std::vector<std::vector<double>> volume_fractions (in.n_evaluation_points(), std::vector<double> (equation_of_state.number_of_lookups()));

      // We need to make a copy of the material model inputs because we want to use the adiabatic pressure
      // rather than the real pressure for the equations of state (to avoid numerical instabilities).
      MaterialModel::MaterialModelInputs<dim> eos_in(in);
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        eos_in.pressure[i] = this->get_adiabatic_conditions().pressure(in.position[i]);

      // Evaluate the equation of state properties over all evaluation points
      equation_of_state.evaluate(eos_in, eos_outputs);
      thermal_conductivity->evaluate(in, out);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0;

          // Calculate volume fractions from mass fractions
          // If there is only one lookup table, set the mass and volume fractions to 1
          std::vector<double> mass_fractions;
          if (equation_of_state.number_of_lookups() == 1)
            mass_fractions.push_back(1.0);
          else
            {
              // We only want to compute mass/volume fractions for fields that are chemical compositions.
              mass_fractions = MaterialUtilities::compute_only_composition_fractions(in.composition[i],
                                                                                     this->introspection().chemical_composition_field_indices());

              // The function compute_volumes_from_masses expects as many mass_fractions as densities.
              // But the function compute_composition_fractions always adds another element at the start
              // of the vector that represents the background field. If there is no lookup table for
              // the background field, the mass_fractions vector is too long and we remove this element.
              if (!has_background_field)
                mass_fractions.erase(mass_fractions.begin());
            }

          volume_fractions[i] = MaterialUtilities::compute_volumes_from_masses(mass_fractions,
                                                                               eos_outputs[i].densities,
                                                                               true);
          if (in.requests_property(MaterialProperties::viscosity) || in.requests_property(MaterialProperties::additional_outputs))
            {
              out.viscosities[i] = std::max(
                                     std::min(viscosity(in.temperature[i], in.pressure[i], volume_fractions[i], in.strain_rate[i],
                                                        in.position[i],in.requests_property(MaterialProperties::additional_outputs),out,i),max_eta),
                                     min_eta
                                   );
            }
          MaterialUtilities::fill_averaged_equation_of_state_outputs(eos_outputs[i], mass_fractions, volume_fractions[i], i, out);
          fill_prescribed_outputs(i, volume_fractions[i], in, out);
        }

      // fill additional outputs if they exist
      equation_of_state.fill_additional_outputs(in, volume_fractions, out);
    }



    template <int dim>
    void
    Steinberger<dim>::
    fill_prescribed_outputs(const unsigned int q,
                            const std::vector<double> &,
                            const MaterialModel::MaterialModelInputs<dim> &in,
                            MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // set up variable to interpolate prescribed field outputs onto compositional field
      const std::shared_ptr<PrescribedFieldOutputs<dim>> prescribed_field_out
        = out.template get_additional_output_object<PrescribedFieldOutputs<dim>>();

      if (this->introspection().composition_type_exists(CompositionalFieldDescription::density)
          && prescribed_field_out != nullptr
          && in.requests_property(MaterialProperties::additional_outputs))
        {
          const unsigned int projected_density_index = this->introspection().find_composition_type(CompositionalFieldDescription::density);
          prescribed_field_out->prescribed_field_outputs[q][projected_density_index] = out.densities[q];
        }
    }



    template <int dim>
    void
    Steinberger<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/steinberger/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
          prm.declare_entry ("Use lateral average temperature for viscosity", "true",
                             Patterns::Bool (),
                             "Whether to use the laterally averaged temperature "
                             "instead of the adiabatic temperature as reference for the "
                             "viscosity calculation. This ensures that the laterally averaged "
                             "viscosities remain more or less constant over the model "
                             "runtime. This behavior might or might not be desired.");
          prm.declare_entry ("Number lateral average bands", "10",
                             Patterns::Integer (1),
                             "Number of bands to compute laterally averaged temperature within.");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Maximum lateral viscosity variation", "1e2",
                             Patterns::Double (0.),
                             "The relative cutoff value for lateral viscosity variations "
                             "caused by temperature deviations. The viscosity may vary "
                             "laterally by this factor squared.");
          prm.declare_entry ("Composition viscosity prefactors", "1",
                             Patterns::Anything (),
                             "List of N prefactors that are used to modify the reference viscosity, "
                             "where N is either equal to one or the number of chemical components "
                             "in the simulation. If only one value is given, then all components "
                             "use the same value. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Viscosity averaging scheme", "harmonic",
                             Patterns::Selection("arithmetic|harmonic|geometric|maximum composition"),
                             "Method to average viscosities over multiple compositional fields. "
                             "One of arithmetic, harmonic, geometric or maximum composition.");

          // Thermal conductivity parameters
          ThermalConductivity::Constant<dim>::declare_parameters (prm);
          prm.declare_entry ("Thermal conductivity formulation", "constant",
                             Patterns::Selection("constant|p-T-dependent"),
                             "Which law should be used to compute the thermal conductivity. "
                             "The 'constant' law uses a constant value for the thermal "
                             "conductivity. The 'p-T-dependent' formulation uses equations "
                             "from Stackhouse et al. (2015): First-principles calculations "
                             "of the lattice thermal conductivity of the lower mantle "
                             "(https://doi.org/10.1016/j.epsl.2015.06.050), and Tosi et al. "
                             "(2013): Mantle dynamics with pressure- and temperature-dependent "
                             "thermal expansivity and conductivity "
                             "(https://doi.org/10.1016/j.pepi.2013.02.004) to compute the "
                             "thermal conductivity in dependence of temperature and pressure. "
                             "The thermal conductivity parameter sets can be chosen in such a "
                             "way that either the Stackhouse or the Tosi relations are used. "
                             "The conductivity description can consist of several layers with "
                             "different sets of parameters. Note that the Stackhouse "
                             "parametrization is only valid for the lower mantle (bridgmanite).");
          ThermalConductivity::TosiStackhouse<dim>::declare_parameters (prm);

          // Table lookup parameters
          EquationOfState::ThermodynamicTableLookup<dim>::declare_parameters(prm);

          // Drucker Prager plasticity parameters
          prm.declare_entry ("Use Drucker-Prager rheology", "false",
                             Patterns::Bool(),
                             "This parameter determines whether to apply plastic yielding "
                             "according to a Drucker-Prager rheology after computing the "
                             "default steinberger viscosity (if true) or not (if false). ");
          Rheology::DruckerPrager<dim>::declare_parameters(prm);

          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    Steinberger<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Steinberger model");
        {
          data_directory                  = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          radial_viscosity_file_name      = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name     = prm.get ("Lateral viscosity file name");
          use_lateral_average_temperature = prm.get_bool ("Use lateral average temperature for viscosity");
          n_lateral_slices                = prm.get_integer("Number lateral average bands");
          min_eta                         = prm.get_double ("Minimum viscosity");
          max_eta                         = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation       = prm.get_double ("Maximum lateral viscosity variation");
          viscosity_averaging_scheme      = MaterialUtilities::parse_compositional_averaging_operation ("Viscosity averaging scheme",
                                            prm);

          // Thermal conductivity parameters
          if (prm.get ("Thermal conductivity formulation") == "constant")
            thermal_conductivity = std::make_unique<ThermalConductivity::Constant<dim>>();
          else if (prm.get ("Thermal conductivity formulation") == "p-T-dependent")
            {
              thermal_conductivity = std::make_unique<ThermalConductivity::TosiStackhouse<dim>>();
              if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(thermal_conductivity.get()))
                sim->initialize_simulator (this->get_simulator());
            }
          else
            AssertThrow(false, ExcMessage("Not a valid thermal conductivity formulation"));

          thermal_conductivity->parse_parameters(prm);

          // Parse the table lookup parameters
          equation_of_state.initialize_simulator (this->get_simulator());
          equation_of_state.parse_parameters(prm);

          // Assign background field and do some error checking
          const unsigned int n_chemical_composition_fields = this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::chemical_composition);
          has_background_field = ((equation_of_state.number_of_lookups() == 1) ||
                                  (equation_of_state.number_of_lookups() == n_chemical_composition_fields + 1));
          AssertThrow ((equation_of_state.number_of_lookups() == 1) ||
                       (equation_of_state.number_of_lookups() == n_chemical_composition_fields) ||
                       (equation_of_state.number_of_lookups() == n_chemical_composition_fields + 1),
                       ExcMessage("The Steinberger material model assumes that either there is a single material "
                                  "in the simulation, or that all compositional fields of the type "
                                  "chemical composition correspond to mass fractions of different materials. "
                                  "There must either be one material lookup file, the same number of "
                                  "material lookup files as compositional fields of type chemical composition, "
                                  "or one additional file (if a background field is used). You have "
                                  + Utilities::int_to_string(equation_of_state.number_of_lookups())
                                  + " material data files, but there are "
                                  + Utilities::int_to_string(n_chemical_composition_fields)
                                  + " fields of type chemical composition."));

          // Parse the Composition viscosity prefactors parameter
          std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
          std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();

          // If there is only one lookup, the list of Viscosity prefactors should have length one.
          // The following if statement applies when defined fields of type chemical composition
          // are not intended to be used as volume fractions of distinct materials.
          if (equation_of_state.number_of_lookups() == 1)
            {
              chemical_field_names.clear();
              compositional_field_names.clear();
            }

          if (has_background_field)
            {
              chemical_field_names.insert(chemical_field_names.begin(),"background");
              compositional_field_names.insert(compositional_field_names.begin(),"background");
            }

          Utilities::MapParsing::Options options(chemical_field_names, "Composition viscosity prefactors");
          options.list_of_allowed_keys = compositional_field_names;
          viscosity_prefactors = Utilities::MapParsing::parse_map_to_double_array(prm.get("Composition viscosity prefactors"), options);

          // Drucker-Prager Plasticity parameters
          enable_drucker_prager_rheology = prm.get_bool ("Use Drucker-Prager rheology");
          drucker_prager_plasticity.initialize_simulator (this->get_simulator());

          drucker_prager_plasticity.parse_parameters(prm);

          prm.leave_subsection();
        }
        prm.leave_subsection();

        // Declare dependencies on solution variables
        this->model_dependence.viscosity = NonlinearDependence::temperature;
        this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.compressibility = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.specific_heat = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
        this->model_dependence.thermal_conductivity = NonlinearDependence::none;
      }
    }



    template <int dim>
    void
    Steinberger<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      equation_of_state.create_additional_named_outputs(out);

      if (this->introspection().composition_type_exists(CompositionalFieldDescription::density)
          && out.template has_additional_output_object<PrescribedFieldOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
      if (enable_drucker_prager_rheology && out.template has_additional_output_object<PlasticAdditionalOutputs<dim>>() == false)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<PlasticAdditionalOutputs<dim>> (n_points));
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Steinberger,
                                   "Steinberger",
                                   "This material model looks up the viscosity from the tables that "
                                   "correspond to the paper of Steinberger and Calderwood "
                                   "2006 (``Models of large-scale viscous flow in the Earth's "
                                   "mantle with constraints from mineral physics and surface observations'', "
                                   "Geophys. J. Int., 167, 1461-1481, "
                                   "<http://dx.doi.org/10.1111/j.1365-246X.2006.03131.x>) and material "
                                   "data from a database generated by the thermodynamics code \\texttt{Perplex}, "
                                   "see <http://www.perplex.ethz.ch/>. "
                                   "The default example data builds upon the thermodynamic "
                                   "database by Stixrude 2011 and assumes a pyrolitic composition by "
                                   "Ringwood 1988 but is easily replaceable by other data files. "
                                   "If 'Use Drucker-Prager rheology' is set to true, this model will "
                                   "apply Drucker-Prager plasticity after computing the default steinberger "
                                   "viscosity. The minimum viscosity with and without the Drucker-Prager "
                                   "plasticity is always equal to the 'Minimum viscosity' in the prm file. ")
  }
}
