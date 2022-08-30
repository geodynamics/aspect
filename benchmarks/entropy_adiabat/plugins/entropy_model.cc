/*
  Copyright (C) 2016 - 2022 by the authors of the ASPECT code.

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


#include "entropy_model.h"
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/table.h>
#include <fstream>
#include <iostream>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    EntropyModel<dim>::initialize()
    {
      AssertThrow (this->get_parameters().formulation_mass_conservation ==
                   Parameters<dim>::Formulation::MassConservation::projected_density_field,
                   ExcMessage("The 'entropy model' material model was only tested with the "
                              "'projected density field' approximation "
                              "for the mass conservation equation, which is not selected."));

      AssertThrow (this->introspection().compositional_name_exists("entropy"),
                   ExcMessage("The 'entropy model' material model requires the existence of a compositional field "
                              "named 'entropy'. This field does not exist."));

      material_lookup = std::make_unique<Utilities::StructuredDataLookup<2>>(7,1.0);
      material_lookup->load_file(data_directory+material_file_name,
                                 this->get_mpi_communicator());
    }



    template <int dim>
    bool
    EntropyModel<dim>::
    is_compressible () const
    {
      return true;
    }



    template <int dim>
    void
    EntropyModel<dim>::evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                                MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      const unsigned int projected_density_index = this->introspection().compositional_index_for_name("density_field");
      const unsigned int entropy_index = this->introspection().compositional_index_for_name("entropy");

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // Use the adiabatic pressure instead of the real one,
          // to stabilize against pressure oscillations in phase transitions.
          // This is a requirement of the projected density approximation for
          // the Stokes equation and not related to the entropy formulation.
          // Also convert pressure from Pa to bar, bar is used in the table.
          Point<2> entropy_pressure(in.composition[i][entropy_index],
                                    this->get_adiabatic_conditions().pressure(in.position[i]) / 1.e5);

          out.densities[i]                      = material_lookup->get_data(entropy_pressure,1);
          out.thermal_expansion_coefficients[i] = material_lookup->get_data(entropy_pressure,2);
          out.specific_heat[i]                  = material_lookup->get_data(entropy_pressure,3);

          const Tensor<1,2> density_gradient    = material_lookup->get_gradients(entropy_pressure,1);
          const Tensor<1,2> pressure_unit_vector ({0.0,1.0});
          out.compressibilities[i]              = (density_gradient * pressure_unit_vector) / out.densities[i];

          // Constant thermal conductivity
          out.thermal_conductivities[i]         = thermal_conductivity_value;

          out.entropy_derivative_pressure[i]    = 0.;
          out.entropy_derivative_temperature[i] = 0.;
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = 0.;

          // set up variable to interpolate prescribed field outputs onto compositional fields
          if (PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>())
            {
              prescribed_field_out->prescribed_field_outputs[i][projected_density_index] = out.densities[i];
            }

          // set up variable to interpolate prescribed field outputs onto temperature field
          if (PrescribedTemperatureOutputs<dim> *prescribed_temperature_out = out.template get_additional_output<PrescribedTemperatureOutputs<dim>>())
            {
              prescribed_temperature_out->prescribed_temperature_outputs[i] = material_lookup->get_data(entropy_pressure,0);
            }

          // Calculate Viscosity
          {
            // read in the viscosity profile
            const double depth = this->get_geometry_model().depth(in.position[i]);
            const double viscosity_profile = depth_dependent_rheology->compute_viscosity(depth);

            // lateral viscosity variations
            const double reference_temperature = this->get_adiabatic_conditions().is_initialized()
                                                 ?
                                                 this->get_adiabatic_conditions().temperature(in.position[i])
                                                 :
                                                 this->get_parameters().adiabatic_surface_temperature;
            const double temperature = material_lookup->get_data(entropy_pressure,0);
            const double delta_temperature = temperature-reference_temperature;

            // Steinberger & Calderwood viscosity
            if (temperature*reference_temperature == 0)
              out.viscosities[i] = min_eta;
            else
              {
                double vis_lateral = std::exp(-lateral_viscosity_prefactor*delta_temperature/(temperature*reference_temperature));
                if (std::isnan(vis_lateral))
                  vis_lateral = 1.0;

                out.viscosities[i] = std::max(std::min(vis_lateral * viscosity_profile,max_eta),min_eta);
              }
          }

          // fill seismic velocities outputs if they exist
          if (SeismicAdditionalOutputs<dim> *seismic_out = out.template get_additional_output<SeismicAdditionalOutputs<dim>>())
            {
              seismic_out->vp[i] = material_lookup->get_data(entropy_pressure,4);
              seismic_out->vs[i] = material_lookup->get_data(entropy_pressure,5);
            }
        }
    }



    template <int dim>
    void
    EntropyModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/opxtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT. ");
          prm.declare_entry ("Material file name", "material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the material data.");
          prm.declare_entry ("Lateral viscosity prefactor", "0.0",
                             Patterns::Double(0),
                             "The lateral viscosity prefactor used for computing the temperature"
                             "dependence of viscosity. This corresponds to H/nR (where H is the "
                             "activation enthalpy, n is the stress exponent and R is the gas "
                             "constant) in an Arrhenius-type viscosity law, and is modeled after "
                             "the law for lateral viscosity variations given in Steinberger and "
                             "Calderwood, 2006. "
                             "\n\n"
                             "Units: $1/K$");
          prm.declare_entry ("Minimum viscosity", "1e19",
                             Patterns::Double (0.),
                             "The minimum viscosity that is allowed in the viscosity "
                             "calculation. Smaller values will be cut off.");
          prm.declare_entry ("Maximum viscosity", "1e23",
                             Patterns::Double (0.),
                             "The maximum viscosity that is allowed in the viscosity "
                             "calculation. Larger values will be cut off.");
          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.leave_subsection();
        }

        // Depth-dependent parameters from the rheology plugin
        Rheology::AsciiDepthProfile<dim>::declare_parameters(prm,
                                                             "Depth dependent viscosity");

        prm.leave_subsection();
      }
    }



    template <int dim>
    void
    EntropyModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Entropy model");
        {
          data_directory              = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_name          = prm.get ("Material file name");
          lateral_viscosity_prefactor = prm.get_double ("Lateral viscosity prefactor");
          min_eta                     = prm.get_double ("Minimum viscosity");
          max_eta                     = prm.get_double ("Maximum viscosity");
          thermal_conductivity_value  = prm.get_double ("Thermal conductivity");

          prm.leave_subsection();
        }

        depth_dependent_rheology = std::make_unique<Rheology::AsciiDepthProfile<dim>>();
        depth_dependent_rheology->initialize_simulator (this->get_simulator());
        depth_dependent_rheology->parse_parameters(prm, "Depth dependent viscosity");
        depth_dependent_rheology->initialize();

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
    EntropyModel<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<SeismicAdditionalOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::SeismicAdditionalOutputs<dim>> (n_points));
        }

      if (out.template get_additional_output<PrescribedFieldOutputs<dim>>() == NULL)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>>
            (n_points, this->n_compositional_fields()));
        }

      if (out.template get_additional_output<PrescribedTemperatureOutputs<dim>>() == NULL)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::PrescribedTemperatureOutputs<dim>>
            (n_points));
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(EntropyModel,
                                   "entropy model",
                                   "A material model that is designed to use pressure and entropy (rather "
                                   "than pressure and temperature) as independent variables. "
                                   "It requires a thermodynamic data table that contains "
                                   "all relevant properties in a specific format as illustrated in "
                                   "the data/material-model/entropy-table/opxtable example folder. "
                                   "The material model requires the use of the projected density "
                                   "approximation for compressibility, and the existence of a "
                                   "compositional field called 'entropy'.")
  }
}
