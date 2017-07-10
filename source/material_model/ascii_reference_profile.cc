/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.

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


#include <aspect/material_model/ascii_reference_profile.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    AsciiReferenceProfile<dim>::AsciiReferenceProfile()
      :
      density_index(numbers::invalid_unsigned_int),
      thermal_expansivity_index(numbers::invalid_unsigned_int),
      specific_heat_index(numbers::invalid_unsigned_int),
      compressibility_index(numbers::invalid_unsigned_int),
      seismic_vp_index(numbers::invalid_unsigned_int),
      seismic_vs_index(numbers::invalid_unsigned_int),
      seismic_dvp_dT_index(numbers::invalid_unsigned_int),
      seismic_dvs_dT_index(numbers::invalid_unsigned_int)
    {}

    template <int dim>
    void
    AsciiReferenceProfile<dim>::initialize ()
    {
      profile.initialize(this->get_mpi_communicator());

      density_index = profile.get_column_index_from_name("density");
      thermal_expansivity_index = profile.get_column_index_from_name("thermal_expansivity");
      specific_heat_index = profile.get_column_index_from_name("specific_heat");
      compressibility_index = profile.get_column_index_from_name("compressibility");

      // these are only optional entries in the data file, read them in if they exist,
      // but keep the invalid unsigned int entry if the columns do not exist
      // the strings have to be all lower case for the lookup class to find them
      seismic_vp_index = profile.maybe_get_column_index_from_name("seismic_vp");
      seismic_vs_index = profile.maybe_get_column_index_from_name("seismic_vs");
      seismic_dvp_dT_index = profile.maybe_get_column_index_from_name("seismic_dvp_dt");
      seismic_dvs_dT_index = profile.maybe_get_column_index_from_name("seismic_dvs_dt");
    }

    template <int dim>
    void
    AsciiReferenceProfile<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature_deviation = in.temperature[i] - this->get_adiabatic_conditions().temperature(position);
          const double pressure_deviation = in.pressure[i] - this->get_adiabatic_conditions().pressure(position);

          const double depth = this->get_geometry_model().depth(position);
          const Point<1> profile_position(depth);

          double visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*temperature_deviation/this->get_adiabatic_conditions().temperature(position)),1e3),1e-3);
          if (std::isnan(visc_temperature_dependence))
            visc_temperature_dependence = 1.0;

          double visc_depth_dependence = viscosity_prefactors[0];
          for (unsigned int j=0; j < transition_depths.size(); ++j)
            {
              if (depth>transition_depths[j])
                visc_depth_dependence = viscosity_prefactors[j+1];
            }

          out.viscosities[i] = viscosity * visc_temperature_dependence * visc_depth_dependence;

          out.thermal_conductivities[i] = thermal_conductivity;

          out.thermal_expansion_coefficients[i] = profile.get_data_component(profile_position,thermal_expansivity_index);
          out.specific_heat[i] = profile.get_data_component(profile_position,specific_heat_index);
          out.compressibilities[i] = profile.get_data_component(profile_position,compressibility_index);

          out.densities[i] = profile.get_data_component(profile_position,density_index)
                             * (1.0 - out.thermal_expansion_coefficients[i] * temperature_deviation)
                             * (tala ? 1.0 : (1.0 + out.compressibilities[i] * pressure_deviation));

          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;

          // fill seismic velocities outputs if they exist
          if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim> >())
            {
              if (seismic_vp_index != numbers::invalid_unsigned_int)
                seismic_out->vp[i] = profile.get_data_component(profile_position,seismic_vp_index);
              if (seismic_vs_index != numbers::invalid_unsigned_int)
                seismic_out->vs[i] = profile.get_data_component(profile_position,seismic_vs_index);
              if (seismic_dvp_dT_index != numbers::invalid_unsigned_int)
                seismic_out->vp[i] += profile.get_data_component(profile_position,seismic_dvp_dT_index)
                                      * temperature_deviation;
              if (seismic_dvs_dT_index != numbers::invalid_unsigned_int)
                seismic_out->vs[i] += profile.get_data_component(profile_position,seismic_dvs_dT_index)
                                      * temperature_deviation;
            }
        }
    }


    template <int dim>
    double
    AsciiReferenceProfile<dim>::
    reference_viscosity () const
    {
      return viscosity;
    }


    template <int dim>
    bool
    AsciiReferenceProfile<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    void
    AsciiReferenceProfile<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Ascii reference profile");
        {
          prm.declare_entry ("Thermal conductivity", "4.0",
                             Patterns::Double (0),
                             "Reference conductivity");
          prm.declare_entry ("Viscosity", "1e21",
                             Patterns::Double (0),
                             "Viscosity");
          prm.declare_entry ("Use TALA", "false",
                             Patterns::Bool (),
                             "Whether to use the TALA instead of the ALA "
                             "approximation.");
          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0),
                             "The temperature dependence of viscosity. Dimensionless exponent.");
          prm.declare_entry ("Transition depths", "1.5e5, 4.1e5, 6.6e5",
                             Patterns::List (Patterns::Double(0)),
                             "A list of depths where the viscosity changes. Values must "
                             "monotonically increase. "
                             "Units: $m$.");
          prm.declare_entry ("Viscosity prefactors", "10, 0.1, 1, 10",
                             Patterns::List (Patterns::Double(0)),
                             "A list of prefactors for the viscosity that determine the viscosity "
                             "profile. Each prefactor is applied in a depth range specified by the "
                             "list of `Transition depths', i.e. the first prefactor is applied above "
                             "the first transition depth, the second one between the first and second "
                             "transition depth, and so on. "
                             "To compute the viscosity profile, this prefactor is multiplied by the "
                             "reference viscosity specified through the parameter `Viscosity'. "
                             "List must have one more entry than Transition depths. "
                             "Units: non-dimensional.");

          aspect::Utilities::AsciiDataProfile<dim>::declare_parameters(prm,
                                                                       "$ASPECT_SOURCE_DIR/data/adiabatic-conditions/ascii-data/",
                                                                       "");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    AsciiReferenceProfile<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Ascii reference profile");
        {
          tala                 = prm.get_bool ("Use TALA");
          thermal_conductivity = prm.get_double ("Thermal conductivity");
          viscosity            = prm.get_double ("Viscosity");
          thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
          transition_depths    = Utilities::string_to_double
                                 (Utilities::split_string_list(prm.get ("Transition depths")));
          viscosity_prefactors = Utilities::string_to_double
                                 (Utilities::split_string_list(prm.get ("Viscosity prefactors")));

          // make sure to check against the depth lists for size errors, since using depth
          if (viscosity_prefactors.size() != transition_depths.size()+1)
            AssertThrow(false, ExcMessage("Error: The list of Viscosity prefactors needs to have exactly "
                                          "one more entry than the list of Transition depths. "));

          profile.parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

    }



    template <int dim>
    void
    AsciiReferenceProfile<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<SeismicAdditionalOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std_cxx11::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
            (new MaterialModel::SeismicAdditionalOutputs<dim> (n_points)));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(AsciiReferenceProfile,
                                   "ascii reference profile",
                                   "A material model that reads in a reference "
                                   "state for density, thermal expansivity, compressibility "
                                   "and specific heat from a text file. "
                                   "\n"
                                   "Note the required format of the "
                                   "input data: The first lines may contain any number of comments "
                                   "if they begin with '#', but one of these lines needs to "
                                   "contain the number of points in the reference state as "
                                   "for example '# POINTS: 3'. "
                                   "Following the comment lines there has to be a single line "
                                   "containing the names of all data columns, separated by arbitrarily "
                                   "many spaces. Column names are not allowed to contain spaces. "
                                   "The file can contain unnecessary columns, but for this plugin it "
                                   "needs to at least provide the columns named `density', "
                                   "`thermal\\_expansivity', `specific\\_heat', and `compressibility'. "
                                   "Note that the data lines in the file need to be sorted in order "
                                   "of increasing depth from 0 to the maximal depth in the model "
                                   "domain. Points in the model that are outside of the provided "
                                   "depth range will be assigned the maximum or minimum depth values, "
                                   "respectively. Points do not need to be equidistant, "
                                   "but the computation of properties is optimized in speed "
                                   "if they are."
                                   "\n"
                                   "\n"
                                   "The viscosity $\\eta$ is computed as "
                                   "\\begin{equation}"
                                   "\\eta(z,T) = \\eta_r(z) \\eta_0 \\exp\\left(-A \\frac{T - T_\\text{adi}}{T_\\text{adi}}\\right),"
                                   "\\end{equation}"
                                   "where $\\eta_r(z)$ is the depth-dependence, which is a "
                                   "piecewise constant function computed according to the "
                                   "list of ``Viscosity prefactors'' and ``Transition depths'', "
                                   "$\\eta_0$ is the reference viscosity specified by the parameter ``Viscosity'' "
                                   "and $A$ describes the dependence on temperature and corresponds to "
                                   "the parameter ``Thermal viscosity exponent''.")
  }
}
