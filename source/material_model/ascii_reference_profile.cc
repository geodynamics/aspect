/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#include <aspect/material_model/ascii_reference_profile.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    AsciiReferenceProfile<dim>::initialize ()
    {
      profile.reset(new aspect::Utilities::AsciiDataProfile<dim>);
      profile->initialize(6,this->get_mpi_communicator());
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
          //const double pressure_deviation = in.pressure[i] - this->get_adiabatic_conditions().pressure(position);

          const double depth = this->get_geometry_model().depth(position);
          const Point<1> profile_position(depth);

          out.viscosities[i] = viscosity;
          out.thermal_conductivities[i] = thermal_conductivity;

          out.specific_heat[i] = profile->get_data_component(profile_position,5);
          out.thermal_expansion_coefficients[i] = profile->get_data_component(profile_position,4);
          out.densities[i] = profile->get_data_component(profile_position,2)
                             * (1.0 - out.thermal_expansion_coefficients[i] * temperature_deviation);

          out.compressibilities[i] = 0.0;
          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
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
    double
    AsciiReferenceProfile<dim>::
    reference_density () const
    {
      return 3340;
    }

    template <int dim>
    double
    AsciiReferenceProfile<dim>::
    reference_cp () const
    {
      return 1200;
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
        }
        prm.leave_subsection();

        aspect::Utilities::AsciiDataProfile<dim>::declare_parameters(prm,
            "$ASPECT_SOURCE_DIR/data/adiabatic-conditions/ascii-data/",
            "simple_test.txt");
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
        }
        prm.leave_subsection();

        profile->parse_parameters(prm);
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::pressure | NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;

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
                                   "A material model for that reads in a reference "
                                   "state from a text file.")
  }
}
