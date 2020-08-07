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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/initial_temperature/continental_geotherm.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/heating_model/interface.h>

#include <cmath>

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    ContinentalGeotherm<dim>::ContinentalGeotherm ()
    {}


    template <int dim>
    void
    ContinentalGeotherm<dim>::
    initialize ()
    {
      // Check that the required radioactive heating model ("compositional heating") is used
      const std::vector<std::string> &heating_models = this->get_heating_model_manager().get_active_heating_model_names();
      AssertThrow(std::find(heating_models.begin(), heating_models.end(), "compositional heating") != heating_models.end(),
                  ExcMessage("The continental geotherm initial temperature plugin requires the compositional heating plugin."));

      // Check that the required material model ("visco plastic") is used
      AssertThrow(Plugins::plugin_type_matches<MaterialModel::ViscoPlastic<dim> >(this->get_material_model()),
                  ExcMessage("The continental geotherm initial temperature plugin requires the viscoplastic material model plugin."));
    }


    template <int dim>
    double
    ContinentalGeotherm<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // Get the depth of the point with respect to the reference surface.
      const double depth = this->get_geometry_model().depth(position);

      // Compute some constants to calculate the temperatures T1 and T2 at the interfaces
      // between the layers upper crust/lower crust (depth y=h1) and lower crust/lithospheric mantle (depth y=h1+h2).
      // T1 = (ab + k2/h2 cdb) / (1 - k2^2/h2^2 db)
      // T2 = (c + k2/h2 T1) d
      // These derivations are based on the assumption that at the boundary between
      // the layers, the heat flows of the two layers are equal, so
      // at y=h1, q1(h1)=q2(h2)=k1 dT1/dy|h1 = k2 dT2/dy|h2 etc.
      const double a = 0.5*densities[0]*heat_productivities[0]*thicknesses[0] + 0.5*densities[1]*heat_productivities[1]*thicknesses[1] + conductivities[0]/thicknesses[0]*T0;
      const double b = 1./(conductivities[0]/thicknesses[0]+conductivities[1]/thicknesses[1]);
      const double c = 0.5*densities[1]*heat_productivities[1]*thicknesses[1] + conductivities[2]/thicknesses[2]*LAB_isotherm;
      const double d = 1./(conductivities[1]/thicknesses[1]+conductivities[2]/thicknesses[2]);

      // Temperature at boundary between layer 1 and 2
      const double T1 = (a*b + conductivities[1]/thicknesses[1]*c*d*b) / (1.-(conductivities[1]*conductivities[1])/(thicknesses[1]*thicknesses[1])*d*b);
      // Temperature at boundary between layer 2 and 3
      const double T2 = (c + conductivities[1]/thicknesses[1]*T1) * d;

      // Temperature in layer 1
      if (depth <= thicknesses[0])
        return -0.5*densities[0]*heat_productivities[0]/conductivities[0]*std::pow(depth,2) + (0.5*densities[0]*heat_productivities[0]*thicknesses[0]/conductivities[0] + (T1-T0)/thicknesses[0])*depth + T0;
      // Temperature in layer 2
      else if (depth <= thicknesses[0]+thicknesses[1])
        return -0.5*densities[1]*heat_productivities[1]/conductivities[1]*std::pow(depth-thicknesses[0],2.) + (0.5*densities[1]*heat_productivities[1]*thicknesses[1]/conductivities[1] + (T2-T1)/thicknesses[1])*(depth-thicknesses[0]) + T1;
      // Temperature in layer 3
      else if (depth <= thicknesses[0]+thicknesses[1]+thicknesses[2])
        return (LAB_isotherm-T2)/thicknesses[2] *(depth-thicknesses[0]-thicknesses[1]) + T2;
      // Return an unrealistically high temperature if we're below the lithosphere.
      // This way we can combine the continental geotherm with an adiabatic profile from the input file
      // using the "minimum" operator on the "List of initial temperature models"
      else
        return std::numeric_limits<double>::max();
    }



    template <int dim>
    void
    ContinentalGeotherm<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Continental geotherm");
        {
          prm.declare_entry ("Layer thicknesses", "30000.",
                             Patterns::List(Patterns::Double(0.)),
                             "List of the 3 thicknesses of the lithospheric layers "
                             "'upper\\_crust', 'lower\\_crust' and 'mantle\\_lithosphere'. "
                             "If only one thickness is given, then the same thickness is used "
                             "for all layers. Units: \\si{meter}.");
          prm.declare_entry ("Surface temperature", "273.15",
                             Patterns::Double (0.),
                             "The value of the surface temperature. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Lithosphere-Asthenosphere boundary isotherm", "1673.15",
                             Patterns::Double (0.),
                             "The value of the isotherm that is assumed at the Lithosphere-"
                             "Asthenosphere boundary. Units: \\si{\\kelvin}.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ContinentalGeotherm<dim>::parse_parameters (ParameterHandler &prm)
    {
      unsigned int n_fields = 0;
      prm.enter_subsection ("Compositional fields");
      {
        n_fields = prm.get_integer ("Number of fields");
      }
      prm.leave_subsection();


      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Continental geotherm");
        {
          LAB_isotherm = prm.get_double ("Lithosphere-Asthenosphere boundary isotherm");
          T0 = prm.get_double ("Surface temperature");
          thicknesses = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Layer thicknesses"))),
                                                                3,
                                                                "Layer thicknesses");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      // Retrieve the indices of the fields that represent the lithospheric layers.
      AssertThrow(this->introspection().compositional_name_exists("upper_crust"),ExcMessage("We need a compositional field called 'upper_crust' representing the upper crust."));
      AssertThrow(this->introspection().compositional_name_exists("lower_crust"),ExcMessage("We need a compositional field called 'lower_crust' representing the lower crust."));
      AssertThrow(this->introspection().compositional_name_exists("lithospheric_mantle"),ExcMessage("We need a compositional field called 'lithospheric_mantle' representing the lithospheric part of the mantle."));

      // For now, we assume a 3-layer system with an upper crust, lower crust and lithospheric mantle
      const unsigned int id_upper_crust = this->introspection().compositional_index_for_name("upper_crust");
      const unsigned int id_lower_crust = this->introspection().compositional_index_for_name("lower_crust");
      const unsigned int id_lithospheric_mantle = this->introspection().compositional_index_for_name("lithospheric_mantle");

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
          heat_productivities.push_back(temp_heat_productivities[id_upper_crust+1]);
          heat_productivities.push_back(temp_heat_productivities[id_lower_crust+1]);
          heat_productivities.push_back(temp_heat_productivities[id_lithospheric_mantle+1]);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Visco Plastic");
        {
          // The material model viscoplastic prefixes an entry for the background material, hence n_fields+1
          const std::vector<double> temp_densities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Densities"))),
                                                     n_fields+1,
                                                     "Densities");
          const std::vector<double> temp_thermal_diffusivities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Thermal diffusivities"))),
                                                                 n_fields+1,
                                                                 "Thermal diffusivities");
          const std::vector<double> temp_heat_capacities = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Heat capacities"))),
                                                           n_fields+1,
                                                           "Heat capacities");

          // The material model viscoplastic prefixes an entry for the background material, hence id+1
          densities.push_back(temp_densities[id_upper_crust+1]);
          densities.push_back(temp_densities[id_lower_crust+1]);
          densities.push_back(temp_densities[id_lithospheric_mantle+1]);

          // Thermal diffusivity kappa = k/(rho*cp), so thermal conductivity k = kappa*rho*cp.
          // The densities are already in the right order, so we don't need to use the compositional
          // field ids.
          conductivities.push_back(temp_thermal_diffusivities[id_upper_crust+1] * densities[0] * temp_heat_capacities[id_upper_crust+1]);
          conductivities.push_back(temp_thermal_diffusivities[id_lower_crust+1] * densities[1] * temp_heat_capacities[id_lower_crust+1]);
          conductivities.push_back(temp_thermal_diffusivities[id_lithospheric_mantle+1] * densities[2] * temp_heat_capacities[id_lithospheric_mantle+1]);

          // To obtain the radioactive heating rate in W/kg, we divide the volumetric heating rate by density
          AssertThrow(heat_productivities.size() == 3 && densities.size() == 3 && conductivities.size() == 3,
                      ExcMessage("The entries for density, conductivity and heat production do not match with the expected number of layers (3)."))

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
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(ContinentalGeotherm,
                                              "continental geotherm",
                                              "This is a temperature initial condition that "
                                              "computes a continental geotherm based on the solution of the "
                                              "steady-state conductive equation $k\\frac{d^2 T}{dy^2}+\\rho H = 0$ "
                                              "as described in e.g. Turcotte and Schubert, "
                                              "Ch. 4.6, or Chapman (1986). As boundary conditions, we take the surface temperature and "
                                              "the temperature of the Lithosphere-Asthenosphere Boundary (LAB). "
                                              "\n"
                                              "The geotherm is computed for a homogeneous lithosphere "
                                              "composed of an upper crust, lower crust and mantle layer. "
                                              "The crustal layers are assumed to have a constant radioactive heating, "
                                              "and all layers are assumed to have a constant thermal conductivity. "
                                              "Layer thicknesses, surface temperature and LAB temperature "
                                              "should be specified by the user. "
                                              "For consistency, the density, heat production and thermal "
                                              "conductivity of each layer are read from the visco plastic material model "
                                              "and the compositional heating model. "
                                              "\n"
                                              "For any depths below the depth of the LAB, a unrealistically high "
                                              "temperature is returned, such that this plugin can be combined with "
                                              "another temperature plugin through the 'minimum' operator. "
                                              "\n"
                                              "Note that the current implementation only works for a 3-layer lithosphere, "
                                              "even though in principle the heat conduction equation can be solved "
                                              "for any number of layers. The naming of the compositional fields "
                                              "that represent the layers is also very specific, namely `upper\\_crust', "
                                              "`lower\\_crust', and `lithospheric\\_mantle'. "
                                              "\n"
                                              "Make sure the top and bottom temperatures of the lithosphere "
                                              "agree with temperatures set in for example the temperature "
                                              "boundary conditions.")
  }
}
