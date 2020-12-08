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


#include <aspect/material_model/steinberger.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>
#include <aspect/lateral_averaging.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>
#include <memory>


namespace aspect
{
  namespace MaterialModel
  {
    namespace internal
    {
      LateralViscosityLookup::LateralViscosityLookup(const std::string &filename,
                                                     const MPI_Comm &comm)
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
            in >> visc;;
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
                                                   const MPI_Comm &comm)
      {
        std::string temp;
        // Read data from disk and distribute among processes
        std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

        min_depth=1e20;
        max_depth=-1;

        while (!in.eof())
          {
            double visc, depth;
            in >> visc;;
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
      // This model allows the user to provide several PerpleX P-T lookup files,
      // each of which corresponds to a different material.
      // The thermodynamic properties and stable mineral phases within each material
      // will in general be different, and must be averaged in a reasonable way.

      // In the following, we populate the unique_phase_names vector.
      // We do this because we are choosing in this material model to combine
      // minerals with different compositions but which are the same phase.
      std::set<std::string> set_phase_volume_column_names;
      for (unsigned i = 0; i < material_file_names.size(); i++)
        {
          material_lookup.push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>
                                    (data_directory+material_file_names[i],interpolation,this->get_mpi_communicator()));

          // Resize the unique_phase_indices object
          unique_phase_indices.resize(material_file_names.size(), std::vector<unsigned int>());

          // Here we look up all of the column names and insert
          // unique names into the unique_phase_names vector and
          // filling the unique_phase_indices object.
          std::vector<std::string> phase_volume_column_names = material_lookup[i]->phase_volume_column_names();
          for (const auto &phase_volume_column_name : phase_volume_column_names)
            {
              // iterate over the present unique_phase_names object
              // to find phase_volume_column_names[j].
              std::vector<std::string>::iterator it = std::find(unique_phase_names.begin(),
                                                                unique_phase_names.end(),
                                                                phase_volume_column_name);

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
                unique_phase_names.push_back(phase_volume_column_name);
            }
        }

      lateral_viscosity_lookup
        = std_cxx14::make_unique<internal::LateralViscosityLookup>(data_directory+lateral_viscosity_file_name,
                                                                   this->get_mpi_communicator());
      radial_viscosity_lookup
        = std_cxx14::make_unique<internal::RadialViscosityLookup>(data_directory+radial_viscosity_file_name,
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
               const std::vector<double> &,
               const SymmetricTensor<2,dim> &,
               const Point<dim> &position) const
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
      const double vis_lateral_exp = -1.0*lateral_viscosity_lookup->lateral_viscosity(depth)*delta_temperature/(temperature*adiabatic_temperature);
      // Limit the lateral viscosity variation to a reasonable interval
      const double vis_lateral = std::max(std::min(std::exp(vis_lateral_exp),max_lateral_eta_variation),1/max_lateral_eta_variation);

      const double vis_radial = radial_viscosity_lookup->radial_viscosity(depth);

      return std::max(std::min(vis_lateral * vis_radial,max_eta),min_eta);
    }



    template <int dim>
    double
    Steinberger<dim>::
    reference_viscosity () const
    {
      return reference_eta;
    }



    template <int dim>
    void
    Steinberger<dim>::
    fill_mass_and_volume_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                    std::vector<std::vector<double>> &mass_fractions,
                                    std::vector<std::vector<double>> &volume_fractions) const
    {
      // Resize mass and volume fraction vectors
      mass_fractions.resize(in.n_evaluation_points(), std::vector<double>(material_lookup.size(), 1.));
      volume_fractions.resize(in.n_evaluation_points(), std::vector<double>(material_lookup.size(), 1.));

      if (material_lookup.size() > 1)
        {
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              double summed_volumes = 0.;

              if (has_background)
                {
                  mass_fractions[i][0] = 1.;
                  for (unsigned int j=1; j<material_lookup.size(); ++j)
                    {
                      const double mass_fraction = in.composition[i][first_composition_index+j-1];
                      mass_fractions[i][j] = mass_fraction;
                      mass_fractions[i][0] -= mass_fraction;
                      volume_fractions[i][j] = mass_fraction/material_lookup[j]->density(in.temperature[i],in.pressure[i]);
                      summed_volumes += volume_fractions[i][j];
                    }
                  volume_fractions[i][0] = mass_fractions[i][0]/material_lookup[0]->density(in.temperature[i],in.pressure[i]);
                  summed_volumes += volume_fractions[i][0];

                }
              else
                {
                  for (unsigned int j=0; j<material_lookup.size(); ++j)
                    {
                      const double mass_fraction = in.composition[i][first_composition_index+j];
                      mass_fractions[i][j] = mass_fraction;
                      volume_fractions[i][j] = mass_fraction/material_lookup[j]->density(in.temperature[i],in.pressure[i]);
                      summed_volumes += volume_fractions[i][j];
                    }
                }

              for (unsigned int j=0; j<material_lookup.size(); ++j)
                volume_fractions[i][j] /= summed_volumes;

            }
        }
    }



    template <int dim>
    void
    Steinberger<dim>::
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
    Steinberger<dim>::
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
    bool
    Steinberger<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    std::array<std::pair<double, unsigned int>,2>
    Steinberger<dim>::
    enthalpy_derivatives (const typename Interface<dim>::MaterialModelInputs &in) const
    {
      // We have to take into account here that the p,T spacing of the table of material properties
      // we use might be on a finer grid than our model. Because of that we compute the enthalpy
      // derivatives by using finite differences that average over the whole temperature and
      // pressure range that is used in this cell. This way we should not miss any phase transformation.
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

      // compute the averaged enthalpy derivatives for all temperatures and
      // pressures in this cell. The 1 means we only do one substep for this
      // computation (see documentation of the called function for more
      // information.
      derivative = material_lookup[0]->enthalpy_derivatives(temperatures,
                                                            pressures,
                                                            1);

      return derivative;
    }



    template <int dim>
    void
    Steinberger<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                               MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      std::vector<std::vector<double>> mass_fractions;
      std::vector<std::vector<double>> volume_fractions;
      fill_mass_and_volume_fractions (in, mass_fractions, volume_fractions);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          if (in.requests_property(MaterialProperties::viscosity))
            out.viscosities[i] = viscosity(in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);

          out.thermal_conductivities[i] = thermal_conductivity_value;
          out.entropy_derivative_pressure[i]    = 0;
          out.entropy_derivative_temperature[i] = 0;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0;

          // The following lines take the appropriate averages of the
          // thermodynamic material properties
          std::vector<double> specific_heats(material_lookup.size());
          std::vector<double> densities(material_lookup.size());
          std::vector<double> compressibilities(material_lookup.size());
          std::vector<double> thermal_expansivities(material_lookup.size());

          for (unsigned int j=0; j<material_lookup.size(); ++j)
            {
              densities[j] = material_lookup[j]->density(in.temperature[i],in.pressure[i]);
              compressibilities[j] = material_lookup[j]->dRhodp(in.temperature[i],in.pressure[i])/densities[j];

              if (!latent_heat)
                {
                  thermal_expansivities[j] = material_lookup[j]->thermal_expansivity(in.temperature[i],in.pressure[i]);
                  specific_heats[j] = material_lookup[j]->specific_heat(in.temperature[i],in.pressure[i]);
                }
            }

          // The density and isothermal compressibility are both volume-averaged
          out.densities[i] = MaterialUtilities::average_value(volume_fractions[i], densities, MaterialUtilities::arithmetic);
          out.compressibilities[i] = MaterialUtilities::average_value(volume_fractions[i], compressibilities, MaterialUtilities::arithmetic);

          if (!latent_heat)
            {
              // Specific heat is measured per unit mass, so it is mass averaged.
              // Thermal expansivity is volume averaged.
              out.specific_heat[i] = MaterialUtilities::average_value(mass_fractions[i], specific_heats, MaterialUtilities::arithmetic);
              out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value(volume_fractions[i], thermal_expansivities, MaterialUtilities::arithmetic);
            }
        }

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
      if (latent_heat)
        {
          double average_temperature(0.0);
          double average_density(0.0);
          std::array<std::pair<double, unsigned int>,2> dH;

          for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
            {
              average_temperature += in.temperature[i];
              average_density += out.densities[i];
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
                  out.thermal_expansion_coefficients[i] = (1 - average_density * dH[1].first) / average_temperature;
                  out.specific_heat[i] = dH[0].first;
                }
              else
                {
                  out.thermal_expansion_coefficients[i] = (1 - out.densities[i] * material_lookup[0]->dHdp(in.temperature[i],pressure)) / in.temperature[i];
                  out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],pressure);
                }
            }
        }

      // fill seismic velocity outputs if they exist
      if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim> >())
        fill_seismic_velocities(in, out.densities, volume_fractions, seismic_out);

      // fill phase volume outputs if they exist
      if (NamedAdditionalMaterialOutputs<dim> *phase_volume_fractions_out = out.template get_additional_output<NamedAdditionalMaterialOutputs<dim> >())
        fill_phase_volume_fractions(in, volume_fractions, phase_volume_fractions_out);
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
                             "assumed to represent the mass fractions of different materials "
                             "and the average property is computed as a sum of "
                             "the value of the compositional field times the "
                             "material property of that field.");
          prm.declare_entry ("Radial viscosity file name", "radial-visc.txt",
                             Patterns::Anything (),
                             "The file name of the radial viscosity data. ");
          prm.declare_entry ("Lateral viscosity file name", "temp-viscosity-prefactor.txt",
                             Patterns::Anything (),
                             "The file name of the lateral viscosity data. ");
          prm.declare_entry ("Use lateral average temperature for viscosity", "true",
                             Patterns::Bool (),
                             "Whether to use to use the laterally averaged temperature "
                             "instead of the adiabatic temperature as reference for the "
                             "viscosity calculation. This ensures that the laterally averaged "
                             "viscosities remain more or less constant over the model "
                             "runtime. This behaviour might or might not be desired.");
          prm.declare_entry ("Number lateral average bands", "10",
                             Patterns::Integer (1),
                             "Number of bands to compute laterally averaged temperature within.");
          prm.declare_entry ("Bilinear interpolation", "true",
                             Patterns::Bool (),
                             "Whether to use bilinear interpolation to compute "
                             "material properties (slower but more accurate). ");
          prm.declare_entry ("Latent heat", "false",
                             Patterns::Bool (),
                             "Whether to include latent heat effects in the "
                             "calculation of thermal expansivity and specific heat. "
                             "Following the approach of Nakagawa et al. 2009. ");
          prm.declare_entry ("Reference viscosity", "1e23",
                             Patterns::Double (0.),
                             "The reference viscosity that is used for pressure scaling. "
                             "To understand how pressure scaling works, take a look at "
                             "\\cite{KHB12}. In particular, the value of this parameter "
                             "would not affect the solution computed by \\aspect{} if "
                             "we could do arithmetic exactly; however, computers do "
                             "arithmetic in finite precision, and consequently we need to "
                             "scale quantities in ways so that their magnitudes are "
                             "roughly the same. As explained in \\cite{KHB12}, we scale "
                             "the pressure during some computations (never visible by "
                             "users) by a factor that involves a reference viscosity. This "
                             "parameter describes this reference viscosity."
                             "\n\n"
                             "For problems with a constant viscosity, you will generally want "
                             "to choose the reference viscosity equal to the actual viscosity. "
                             "For problems with a variable viscosity, the reference viscosity "
                             "should be a value that adequately represents the order of "
                             "magnitude of the viscosities that appear, such as an average "
                             "value or the value one would use to compute a Rayleigh number."
                             "\n\n"
                             "Units: \\si{\\pascal\\second}.");
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
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
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
          data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_names  = Utilities::split_string_list
                                 (prm.get ("Material file names"));
          radial_viscosity_file_name   = prm.get ("Radial viscosity file name");
          lateral_viscosity_file_name  = prm.get ("Lateral viscosity file name");
          use_lateral_average_temperature = prm.get_bool ("Use lateral average temperature for viscosity");
          n_lateral_slices = prm.get_integer("Number lateral average bands");
          interpolation        = prm.get_bool ("Bilinear interpolation");
          latent_heat          = prm.get_bool ("Latent heat");
          reference_eta        = prm.get_double ("Reference viscosity");
          min_eta              = prm.get_double ("Minimum viscosity");
          max_eta              = prm.get_double ("Maximum viscosity");
          max_lateral_eta_variation    = prm.get_double ("Maximum lateral viscosity variation");
          thermal_conductivity_value = prm.get_double ("Thermal conductivity");

          prm.leave_subsection();
        }
        prm.leave_subsection();

        // Do some error checking
        AssertThrow ((material_file_names.size() == 1) ||
                     (material_file_names.size() == this->n_compositional_fields()) ||
                     (material_file_names.size() == this->n_compositional_fields() + 1),
                     ExcMessage("This material model expects either one material data file, or as many files as compositional fields, "
                                "or as many files as compositional fields plus one (in which case the first file "
                                "is assumed to contain a background composition). This condition is not fulfilled. You "
                                "prescribed " + Utilities::int_to_string(material_file_names.size()) + " material data files, but there are " +
                                Utilities::int_to_string(this->n_compositional_fields()) + " compositional fields."));

        // The Steinberger material model currently assumes that all the
        // compositional fields correspond to materials with lookup tables.
        // Therefore the first composition index is hard-coded as zero.
        // If more compositional fields are needed, this parameter allows
        // the user to declare the first compositional index which
        // corresponds to a lookup.
        first_composition_index = 0;

        // Define whether there is a background compositional field
        has_background = (material_file_names.size() == this->n_compositional_fields() + 1);

        if (latent_heat)
          AssertThrow (material_file_names.size() == 1,
                       ExcMessage("Isochemical latent heat calculations are only implemented for a single material lookup."));

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
                                   "\\url{http://dx.doi.org/10.1111/j.1365-246X.2006.03131.x}) and material "
                                   "data from a database generated by the thermodynamics code \\texttt{Perplex}, "
                                   "see \\url{http://www.perplex.ethz.ch/}. "
                                   "The default example data builds upon the thermodynamic "
                                   "database by Stixrude 2011 and assumes a pyrolitic composition by "
                                   "Ringwood 1988 but is easily replaceable by other data files. ")
  }
}
