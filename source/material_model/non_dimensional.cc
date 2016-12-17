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


#include <aspect/material_model/non_dimensional.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Nondimensional<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          const Point<dim> position = in.position[i];
          const double temperature_deviation = in.temperature[i] - this->get_adiabatic_conditions().temperature(position);
          const double pressure_deviation = in.pressure[i] - this->get_adiabatic_conditions().pressure(position);

          const double depth = this->get_geometry_model().depth(position);

          out.viscosities[i] = (compressible ? Di : 1.0) / Ra
                               * exp(-b*temperature_deviation + c*depth);

          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = 1.0;
          out.thermal_expansion_coefficients[i] = compressible ? Di : 1.0;

          out.densities[i] = reference_rho
                             * exp(depth * Di/gamma)
                             * (1.0
                                - out.thermal_expansion_coefficients[i] * temperature_deviation
                                + (tala?0.0:1.0)*Di*gamma * pressure_deviation );

          out.compressibilities[i] = 0.0;
          out.entropy_derivative_pressure[i] = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }
    }


    template <int dim>
    double
    Nondimensional<dim>::
    reference_viscosity () const
    {
      return compressible ? (Di/Ra) : (1.0/Ra);
    }

    template <int dim>
    double
    Nondimensional<dim>::
    reference_density () const
    {
      return 1;
    }

    template <int dim>
    double
    Nondimensional<dim>::
    reference_cp () const
    {
      return reference_specific_heat;
    }

    template <int dim>
    bool
    Nondimensional<dim>::
    is_compressible () const
    {
      return compressible;
    }



    template <int dim>
    void
    Nondimensional<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Nondimensional model");
        {
          prm.declare_entry ("Reference density", "1.0",
                             Patterns::Double (0),
                             "Reference density $\\rho_0$. Units: $kg/m^3$.");
          prm.declare_entry ("Ra", "1e4",
                             Patterns::Double (0),
                             "Rayleigh number Ra");
          prm.declare_entry ("Di", "0.0",
                             Patterns::Double (0),
                             "Dissipation number. Pick 0.0 for incompressible "
                             "computations.");
          prm.declare_entry ("gamma", "1.0",
                             Patterns::Double (0),
                             "Grueneisen parameter");
          prm.declare_entry ("Reference specific heat", "1.0",
                             Patterns::Double (0),
                             "The value of the specific heat $cp$. "
                             "Units: $J/kg/K$.");
          prm.declare_entry ("b", "0.0",
                             Patterns::Double (0),
                             "Temperature prefactor for viscosity.");
          prm.declare_entry ("c", "0.0",
                             Patterns::Double (0),
                             "Depth prefactor for viscosity.");
          prm.declare_entry ("Use TALA", "false",
                             Patterns::Bool (),
                             "Whether to use the TALA instead of the ALA "
                             "approximation.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Nondimensional<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Nondimensional model");
        {
          reference_rho   = prm.get_double ("Reference density");
          b               = prm.get_double ("b");
          c               = prm.get_double ("c");
          Di              = prm.get_double ("Di");
          Ra              = prm.get_double ("Ra");
          gamma           = prm.get_double ("gamma");

          tala            = prm.get_bool ("Use TALA");

          reference_specific_heat = prm.get_double ("Reference specific heat");
          compressible = (Di!=0.0);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature;
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
    ASPECT_REGISTER_MATERIAL_MODEL(Nondimensional,
                                   "nondimensional",
                                   "A material model for nondimensionalized "
                                   "computations for compressible or incompressible "
                                   "computations defined through Rayleigh number Ra "
                                   "and Dissipation number Di. This model is made "
                                   "to be used with the Boussinesq, ALA, or TALA "
                                   "formulation.")
  }
}
