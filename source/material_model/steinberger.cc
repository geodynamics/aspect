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
      std::set<std::string> set_phase_names;
      for (unsigned i = 0; i < material_file_names.size(); i++)
        {
          material_lookup.push_back(std_cxx14::make_unique<MaterialModel::MaterialUtilities::Lookup::PerplexReader>
                                    (data_directory+material_file_names[i],interpolation,this->get_mpi_communicator()));

          // Here we lookup all of the different phases and insert them into the set_phase_names set.
          std::vector<std::string> phase_names = material_lookup[i]->phase_volume_phase_names();
          unique_phase_indices.push_back(std::vector<int>());
          for (std::string const &name : phase_names)
            {
              set_phase_names.insert(name);
            }
        }
      // Copy the set of phase names into the unique_phase_names vector
      std::copy(set_phase_names.begin(), set_phase_names.end(), std::back_inserter(unique_phase_names));

      // We now loop over all the material file names again, finding the
      // index of material_lookup[i]->phase_volume_phase_names()
      for (unsigned i = 0; i < material_file_names.size(); i++)
        {
          std::vector<std::string> phase_names = material_lookup[i]->phase_volume_phase_names();
          for (unsigned j = 0; j < unique_phase_names.size(); j++)
            {
              unique_phase_indices[i].push_back(material_lookup[i]->phase_volume_index(unique_phase_names[j]));
            }
        }

      lateral_viscosity_lookup
        = std_cxx14::make_unique<internal::LateralViscosityLookup>(data_directory+lateral_viscosity_file_name,
                                                                   this->get_mpi_communicator());
      radial_viscosity_lookup
        = std_cxx14::make_unique<internal::RadialViscosityLookup>(data_directory+radial_viscosity_file_name,
                                                                  this->get_mpi_communicator());
      avg_temp.resize(n_lateral_slices);

    }



    template <int dim>
    void
    Steinberger<dim>::
    update()
    {
      if (use_lateral_average_temperature)
        {
          this->get_lateral_averaging().get_temperature_averages(avg_temp);
          for (unsigned int i = 0; i < avg_temp.size(); ++i)
            AssertThrow(numbers::is_finite(avg_temp[i]),
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
          const unsigned int idx = static_cast<unsigned int>((avg_temp.size()-1) * depth / this->get_geometry_model().maximal_depth());
          delta_temperature = temperature-avg_temp[idx];
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
    double
    Steinberger<dim>::
    specific_heat (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &) const
    {
      double cp = 0.0;

      if (material_lookup.size() == 1)
        {
          cp = material_lookup[0]->specific_heat(temperature,pressure);
        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          const double background_cp = material_lookup[0]->specific_heat(temperature,pressure);
          cp = background_cp;
          for (unsigned int i = 0; i < compositional_fields.size(); ++i)
            cp += compositional_fields[i] *
                  (material_lookup[i+1]->specific_heat(temperature,pressure) - background_cp);
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); ++i)
            cp += compositional_fields[i] * material_lookup[i]->specific_heat(temperature,pressure);
        }

      return cp;
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &,
                          const Point<dim> &) const
    {
      return thermal_conductivity_value;
    }



    template <int dim>
    double
    Steinberger<dim>::
    density (const double temperature,
             const double pressure,
             const std::vector<double> &compositional_fields,
             const Point<dim> &) const
    {
      double rho = 0.0;
      if (material_lookup.size() == 1)
        {
          rho = material_lookup[0]->density(temperature,pressure);
        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          const double background_density = material_lookup[0]->density(temperature,pressure);
          rho = background_density;
          for (unsigned int i = 0; i < compositional_fields.size(); ++i)
            rho += compositional_fields[i] *
                   (material_lookup[i+1]->density(temperature,pressure) - background_density);
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); ++i)
            rho += compositional_fields[i] * material_lookup[i]->density(temperature,pressure);
        }

      return rho;
    }



    template <int dim>
    double
    Steinberger<dim>::
    thermal_expansion_coefficient (const double      temperature,
                                   const double      pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &/*position*/) const
    {
      double alpha = 0.0;

      if (material_lookup.size() == 1)
        {
          alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          const double background_alpha = material_lookup[0]->thermal_expansivity(temperature,pressure);
          alpha = background_alpha;
          for (unsigned int i = 0; i<compositional_fields.size(); ++i)
            alpha += compositional_fields[i] *
                     (material_lookup[i+1]->thermal_expansivity(temperature,pressure) - background_alpha);
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); ++i)
            alpha += compositional_fields[i] * material_lookup[i]->thermal_expansivity(temperature,pressure);
        }
      return alpha;
    }



    template <int dim>
    double
    Steinberger<dim>::
    seismic_Vp (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &) const
    {
      double vp = 0.0;

      if (material_lookup.size() == 1)
        {
          vp = material_lookup[0]->seismic_Vp(temperature,pressure);
        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          const double background_vp = material_lookup[0]->seismic_Vp(temperature,pressure);
          vp = background_vp;
          for (unsigned int i = 0; i < compositional_fields.size(); ++i)
            vp += compositional_fields[i] *
                  (material_lookup[i+1]->seismic_Vp(temperature,pressure) - background_vp);
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); i++)
            vp += compositional_fields[i] * material_lookup[i]->seismic_Vp(temperature,pressure);
        }
      return vp;
    }



    template <int dim>
    double
    Steinberger<dim>::
    seismic_Vs (const double      temperature,
                const double      pressure,
                const std::vector<double> &compositional_fields,
                const Point<dim> &) const
    {
      double vs = 0.0;

      if (material_lookup.size() == 1)
        {
          vs = material_lookup[0]->seismic_Vs(temperature,pressure);
        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          const double background_vs = material_lookup[0]->seismic_Vs(temperature,pressure);
          vs = background_vs;
          for (unsigned int i = 0; i < compositional_fields.size(); ++i)
            vs += compositional_fields[i] *
                  (material_lookup[i+1]->seismic_Vs(temperature,pressure) - background_vs);
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); i++)
            vs += compositional_fields[i] * material_lookup[i]->seismic_Vs(temperature,pressure);
        }
      return vs;
    }


    template <int dim>
    std::vector<double>
    Steinberger<dim>::
    phase_volume_fractions (const double      temperature,
                            const double      pressure,
                            const std::vector<double> &compositional_fields,
                            const Point<dim> &) const
    {
      std::vector<double> volume_fractions(unique_phase_names.size(), 0.);

      if (material_lookup.size() == 1)
        {
          for (unsigned int j = 0; j < unique_phase_names.size(); j++)
            {
              if (unique_phase_indices[0][j] != -1)
                volume_fractions[j] += material_lookup[0]->phase_volume_fraction(unique_phase_indices[0][j],temperature,pressure);
            }

        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          for (unsigned int j = 0; j < unique_phase_names.size(); j++)
            {
              double background_volume_fraction = 0.;
              if (unique_phase_indices[0][j] != -1)
                background_volume_fraction = material_lookup[0]->phase_volume_fraction(unique_phase_indices[0][j],temperature,pressure);

              volume_fractions[j] = background_volume_fraction;
              for (unsigned int i = 0; i < compositional_fields.size(); ++i)
                {
                  if (unique_phase_indices[i+1][j] != -1)
                    volume_fractions[j] += compositional_fields[i] * (material_lookup[i+1]->phase_volume_fraction(unique_phase_indices[i+1][j],temperature,pressure) - background_volume_fraction);
                }
            }
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); i++)
            {
              for (unsigned int j = 0; j < unique_phase_names.size(); j++)
                {
                  if (unique_phase_indices[i][j] != -1)
                    volume_fractions[j] += compositional_fields[i] * material_lookup[i]->phase_volume_fraction(unique_phase_indices[i][j],temperature,pressure);
                }
            }
        }
      return volume_fractions;
    }



    template <int dim>
    double
    Steinberger<dim>::
    compressibility (const double temperature,
                     const double pressure,
                     const std::vector<double> &compositional_fields,
                     const Point<dim> &position) const
    {
      double dRhodp = 0.0;
      if (material_lookup.size() == 1)
        {
          dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
        }
      else if (material_lookup.size() == compositional_fields.size() + 1)
        {
          const double background_dRhodp = material_lookup[0]->dRhodp(temperature,pressure);
          dRhodp = background_dRhodp;
          for (unsigned int i = 0; i < compositional_fields.size(); ++i)
            dRhodp += compositional_fields[i] *
                      (material_lookup[i+1]->dRhodp(temperature,pressure) - background_dRhodp);
        }
      else
        {
          for (unsigned i = 0; i < material_lookup.size(); i++)
            dRhodp += compositional_fields[i] * material_lookup[i]->dRhodp(temperature,pressure);
        }

      const double rho = density(temperature,pressure,compositional_fields,position);
      return (1/rho)*dRhodp;
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
    enthalpy_derivative (const typename Interface<dim>::MaterialModelInputs &in) const
    {
      // We have to take into account here that the p,T spacing of the table of material properties
      // we use might be on a finer grid than our model. Because of that we compute the enthalpy
      // derivatives by using finite differences that average over the whole temperature and
      // pressure range that is used in this cell. This way we should not miss any phase transformation.
      std::array<std::pair<double, unsigned int>,2> derivative;

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
      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // We are only asked to give viscosities if strain_rate.size() > 0.
          if (in.requests_property(MaterialProperties::viscosity))
            out.viscosities[i]                  = viscosity                     (in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);

          out.densities[i]                      = density                       (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          if (!latent_heat)
            {
              out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
              out.specific_heat[i]                  = specific_heat                 (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
            }
          out.thermal_conductivities[i]         = thermal_conductivity          (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.compressibilities[i]              = compressibility               (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.entropy_derivative_pressure[i]    = 0;
          out.entropy_derivative_temperature[i] = 0;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0;

          // fill seismic velocity outputs if they exist
          if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim> >())
            {
              seismic_out->vp[i] = seismic_Vp(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
              seismic_out->vs[i] = seismic_Vs(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
            }

          // fill phase volume outputs if they exist
          if (PhasePropertyOutputs<dim> *phase_volume_fractions_out = out.template get_additional_output<PhasePropertyOutputs<dim> >())
            {
              phase_volume_fractions_out->phase_property_values[i] = phase_volume_fractions(in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
            }
        }

      if (latent_heat)
        {
          /* We separate the calculation of specific heat and thermal expansivity,
           * because they may depend on cell-wise averaged values that are only
           * available here.
           */
          double average_temperature(0.0);
          double average_density(0.0);
          for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
            {
              average_temperature += in.temperature[i];
              average_density += out.densities[i];
            }
          average_temperature /= in.n_evaluation_points();
          average_density /= in.n_evaluation_points();

          std::array<std::pair<double, unsigned int>,2> dH;
          if (in.current_cell.state() == IteratorState::valid)
            dH = enthalpy_derivative(in);

          for (unsigned int i = 0; i < in.n_evaluation_points(); ++i)
            {
              // Use the adiabatic pressure instead of the real one,
              // to stabilize against pressure oscillations in phase transitions
              const double pressure = this->get_adiabatic_conditions().pressure(in.position[i]);

              // If all of the derivatives were computed successfully
              if ((in.current_cell.state() == IteratorState::valid)
                  && (dH[0].second > 0)
                  && (dH[1].second > 0))
                {
                  // alpha = (1 - rho * dH/dp) / T
                  out.thermal_expansion_coefficients[i] = (1 - average_density * dH[1].first) / average_temperature;
                  // cp = dH/dT
                  out.specific_heat[i] = dH[0].first;
                }
              else
                {
                  out.thermal_expansion_coefficients[i] = (1 - out.densities[i] * material_lookup[0]->dHdp(in.temperature[i],pressure)) / in.temperature[i];
                  out.specific_heat[i] = material_lookup[0]->dHdT(in.temperature[i],pressure);
                }
            }
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

        if (latent_heat)
          AssertThrow (material_file_names.size() == 1,
                       ExcMessage("This formalism is currently only implemented for one material "
                                  "table."));

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
      if (out.template get_additional_output<SeismicAdditionalOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (out.template get_additional_output<PhasePropertyOutputs<dim> >() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::PhasePropertyOutputs<dim>> (n_points, unique_phase_names));
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
