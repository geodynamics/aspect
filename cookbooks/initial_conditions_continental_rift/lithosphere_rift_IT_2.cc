/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include "lithosphere_rift_IT.h"
#include <aspect/utilities.h>
#include "lithosphere_rift_IC.h"
#include <aspect/geometry_model/box.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/heating_model/interface.h>

#include <cmath>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    LithosphereRift<dim>::LithosphereRift ()
    {}


    template <int dim>
    void
    LithosphereRift<dim>::
    initialize ()
    {

      // Check that the required radioactive heating model ("compositional heating") is used.
      const std::vector<std::string> &heating_models = this->get_heating_model_manager().get_active_heating_model_names();
      AssertThrow(std::find(heating_models.begin(), heating_models.end(), "compositional heating") != heating_models.end(),
                  ExcMessage("The lithosphere with rift initial temperature plugin requires the compositional heating plugin."));

      // Check that the required material model ("visco plastic") is used.
      AssertThrow(Plugins::plugin_type_matches<const MaterialModel::ViscoPlastic<dim>> (this->get_material_model()),
                  ExcMessage("The lithosphere with rift initial temperature plugin requires the viscoplastic material model plugin."));

      // The simulator only keeps the initial conditions around for
      // the first time step. As a consequence, we have to save a
      // shared pointer to that object ourselves the first time we get here.
      if (initial_composition_manager == nullptr)
        initial_composition_manager = this->get_initial_composition_manager_pointer();

      // Check that the required initial composition model is used.
      // This doesn't work during initialization.
      //AssertThrow(initial_composition_manager->template has_matching_initial_composition_model<const InitialComposition::LithosphereRift<dim>>(),
      //             ExcMessage("The initial temperature plugin lithosphere with rift requires the corresponding initial composition plugin."));

      // Determine whether a cartesian or a spherical geometry is used.
      cartesian_geometry = Plugins::plugin_type_matches<const GeometryModel::Box<dim>> (this->get_geometry_model());
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      Assert(initial_composition_manager->template has_matching_initial_composition_model<const InitialComposition::LithosphereRift<dim>>(),
             ExcMessage("The initial temperature plugin lithosphere with rift requires the corresponding initial composition plugin."));

      // Get the initial composition plugin
      const InitialComposition::LithosphereRift<dim> &ic = initial_composition_manager->template get_matching_initial_composition_model<const InitialComposition::LithosphereRift<dim>>();

      // Convert to surface position
      const Point<dim - 1> surface_position = ic.surface_position(position, cartesian_geometry);

      // Compute the local thickness of the upper crust, lower crust and mantle part of the lithosphere
      // based on the distance from the rift axis and polygon edges.
      std::vector<double> local_thicknesses = ic.compute_local_thicknesses(surface_position);

      // Get depth with respect to the surface.
      const double depth = this->get_geometry_model().depth(position);

      return temperature(depth, local_thicknesses);
    }


    template <int dim>
    double
    LithosphereRift<dim>::
    temperature (const double depth,
                 const std::vector<double> layer_thicknesses) const
    {
      // Compute some constants
      const double a = 0.5*densities[0]*heat_productivities[0]*layer_thicknesses[0] + 0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1] + conductivities[0]/layer_thicknesses[0]*T0;
      const double b = 1./(conductivities[0]/layer_thicknesses[0]+conductivities[1]/layer_thicknesses[1]);
      const double c = 0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1] + conductivities[2]/layer_thicknesses[2]*LAB_isotherm;
      const double d = 1./(conductivities[1]/layer_thicknesses[1]+conductivities[2]/layer_thicknesses[2]);

      //Temperature at boundary between layer 1 and 2
      const double T1 = (a*b + conductivities[1]/layer_thicknesses[1]*c*d*b) / (1.-(conductivities[1]*conductivities[1])/(layer_thicknesses[1]*layer_thicknesses[1])*d*b);
      //Temperature at boundary between layer 2 and 3
      const double T2 = (c + conductivities[1]/layer_thicknesses[1]*T1) * d;

      // Temperature in layer 1
      if (depth < layer_thicknesses[0])
        return -0.5*densities[0]*heat_productivities[0]/conductivities[0]*std::pow(depth,2) + (0.5*densities[0]*heat_productivities[0]*layer_thicknesses[0]/conductivities[0] + (T1-T0)/layer_thicknesses[0])*depth + T0;
      // Temperature in layer 2
      else if (depth < layer_thicknesses[0]+layer_thicknesses[1])
        return -0.5*densities[1]*heat_productivities[1]/conductivities[1]*std::pow(depth-layer_thicknesses[0],2.) + (0.5*densities[1]*heat_productivities[1]*layer_thicknesses[1]/conductivities[1] + (T2-T1)/layer_thicknesses[1])*(depth-layer_thicknesses[0]) + T1;
      // Temperature in layer 3
      else if (depth < layer_thicknesses[0]+layer_thicknesses[1]+layer_thicknesses[2])
        return (LAB_isotherm-T2)/layer_thicknesses[2] *(depth-layer_thicknesses[0]-layer_thicknesses[1]) + T2;
      // Temperature in the sublithospheric mantle up to the user-set compensation depth
      // equals the temperature at the lithosphere-asthenosphere boundary.
      else if (use_compensation_depth && depth < compensation_depth)
        return LAB_isotherm;
      // Return a NaN sublithospheric temperature.
      // This way we can combine the continental geotherm with an adiabatic profile from the input file
      // using the "if valid" operator on the "List of initial temperature models".
      else
        return std::numeric_limits<double>::quiet_NaN();
    }


    template <int dim>
    void
    LithosphereRift<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          prm.declare_entry("Surface temperature", "273.15",
                            Patterns::Double(0),
                            "The value of the surface temperature. Units: \\si{\\kelvin}.");
          prm.declare_entry("LAB isotherm temperature", "1673.15",
                            Patterns::Double(0),
                            "The value of the isothermal boundary temperature assumed at the LAB "
                            "and up to the compensation depth when 'Use temperature compensation depth' "
                            "is set to true. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Use temperature compensation depth", "false",
                             Patterns::Bool (),
                             "Whether or not to use a compensation depth up to which the LAB isotherm temperature "
                             "is prescribed below the lithosphere. If true, this plugin can be combined with "
                             "the adabiatic temperature plugin. The adiabatic surface temperature can be adapted to "
                             "ensure the LAB temperature and the adiabatic temperature match at the compensation depth. "
                             "If false, combining this plugin with the adiabatic temperature plugin will lead to "
                             "jumps in the temperature profile at the LAB when the lithospheric thickness varies "
                             "laterally. Units: -.");
          prm.declare_entry("Temperature compensation depth", "120e3",
                            Patterns::Double(0),
                            "The depth of the temperature compensation, i.e. the depth up to which the LAB isotherm "
                            "is prescribed below the lithosphere. Units: \\si{\\meter}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    LithosphereRift<dim>::parse_parameters (ParameterHandler &prm)
    {
      unsigned int n_fields = 0;
      prm.enter_subsection ("Compositional fields");
      {
        n_fields = prm.get_integer ("Number of fields");
      }
      prm.leave_subsection();

      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Lithosphere with rift");
        {
          LAB_isotherm = prm.get_double ("LAB isotherm temperature");
          T0 = prm.get_double ("Surface temperature");
          use_compensation_depth = prm.get_bool ("Use temperature compensation depth");
          compensation_depth = prm.get_double ("Temperature compensation depth");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Retrieve the indices of the fields that represent the lithospheric layers.
      AssertThrow(this->introspection().compositional_name_exists("upper"),ExcMessage("We need a compositional field called 'upper' representing the upper crust."));
      AssertThrow(this->introspection().compositional_name_exists("lower"),ExcMessage("We need a compositional field called 'lower' representing the lower crust."));
      AssertThrow(this->introspection().compositional_name_exists("mantle_L"),ExcMessage("We need a compositional field called 'mantle_L' representing the lithospheric part of the mantle."));

      // For now, we assume a 3-layer system with an upper crust, lower crust and lithospheric mantle
      const unsigned int id_upper = this->introspection().compositional_index_for_name("upper");
      const unsigned int id_lower = this->introspection().compositional_index_for_name("lower");
      const unsigned int id_mantle_L = this->introspection().compositional_index_for_name("mantle_L");

      // Retrieve other material properties set in different sections such that there
      // is no need to set them twice.

      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Compositional heating");
        {
          // The heating model compositional heating prefixes an entry for the background material
          const std::vector<double> temp_heat_productivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Compositional heating values"))),
                                                               n_fields+1,
                                                               "Compositional heating values");
          // This sets the heat productivity in W/m3 units
          heat_productivities.push_back(temp_heat_productivities[id_upper+1]);
          heat_productivities.push_back(temp_heat_productivities[id_lower+1]);
          heat_productivities.push_back(temp_heat_productivities[id_mantle_L+1]);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection ("Compositional fields");
      {
        list_of_composition_names = Utilities::split_string_list (prm.get("Names of fields"));
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco Plastic");
        {
          const std::vector<double> temp_thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                 n_fields+1,
                                                                 "Thermal diffusivities");
          const std::vector<double> temp_heat_capacities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacities"))),
                                                           n_fields+1,
                                                           "Heat capacities");

          n_phases_for_each_composition = std::make_unique<std::vector<unsigned int>>();
          const std::vector<double> temp_densities = Utilities::parse_map_to_double_array (prm.get("Densities"),
                                                     list_of_composition_names,
                                                     /*has_background_field=*/true,
                                                     "Densities",
                                                     true,
                                                     n_phases_for_each_composition);

          // Assemble a list of phase densities for each composition.
          // Add 1 for background material.
          std::vector<std::vector<double>> densities_per_composition(this->n_compositional_fields()+1);
          unsigned int counter = 0;
          for (unsigned int i = 0; i < (*n_phases_for_each_composition).size(); ++i)
            {
              for (unsigned int j = 0; j < (*n_phases_for_each_composition)[i]; ++j)
                {
                  densities_per_composition[i].push_back(temp_densities[counter]);
                  ++counter;
                }
            }

          densities.push_back(densities_per_composition[id_upper+1][0]);
          densities.push_back(densities_per_composition[id_lower+1][0]);
          densities.push_back(densities_per_composition[id_mantle_L+1][0]);

          // Thermal diffusivity kappa = k/(rho*cp), so thermal conducitivity k = kappa*rho*cp
          conductivities.push_back(temp_thermal_diffusivities[id_upper+1] * densities[0] * temp_heat_capacities[id_upper+1]);
          conductivities.push_back(temp_thermal_diffusivities[id_lower+1] * densities[1] * temp_heat_capacities[id_lower+1]);
          conductivities.push_back(temp_thermal_diffusivities[id_mantle_L+1] * densities[2] * temp_heat_capacities[id_mantle_L+1]);

          // To obtain the radioactive heating rate in W/kg, we divide the volumetric heating rate by density
          AssertThrow(heat_productivities.size() == 3 && densities.size() == 3 && conductivities.size() == 3,
                      ExcMessage("The entries for density, conductivity and heat production do not match with the expected number of layers (3)."));

          for (unsigned int i = 0; i<3; ++i)
            heat_productivities[i] /= densities[i];
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }


}




// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(LithosphereRift,
                                              "lithosphere with rift",
                                              "This is a temperature initial condition that "
                                              "computes a continental geotherm based on the "
                                              "conductive equations of Turcotte and Schubert Ch. 4.6. "
                                              "This plugin only works with the corresponding composition "
                                              "initial conditions plugin 'lithosphere with rift'. "
                                              "The geotherm can be computed for any number of crustal "
                                              "layers, for each of which a density, heat production and thermal "
                                              "conductivity should be supplied. "
                                              "Make sure the top and bottom temperatures of the lithosphere "
                                              "agree with temperatures set in for example the temperature "
                                              "boundary conditions. "
                                              "The thickness of the lithospheric layers is adapted for distance "
                                              "to the rift specified with the InitialComposition plugin `lithosphere "
                                              "with rift'. Additional variations in lithospheric thickness can be "
                                              "added through specifying polygonal areas of different thicknesses. ")
  }
}
