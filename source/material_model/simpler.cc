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


#include <aspect/material_model/simpler.h>
#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    bool
    Simpler<dim>::
    is_compressible () const
    {
      return equation_of_state.is_compressible ();
    }

    template <int dim>
    double
    Simpler<dim>::
    reference_viscosity () const
    {
      return constant_rheology.compute_viscosity();
    }

    template <int dim>
    void
    Simpler<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // The Simpler model does not depend on composition
      EquationOfStateOutputs<dim> eos_outputs (1);

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs);

          out.viscosities[i] = constant_rheology.compute_viscosity();
          out.densities[i] = eos_outputs.densities[0];
          out.thermal_expansion_coefficients[i] = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[i] = eos_outputs.specific_heat_capacities[0];
          out.thermal_conductivities[i] = k_value;
          out.compressibilities[i] = eos_outputs.compressibilities[0];

          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c] = 0.0;
        }

    }


    template <int dim>
    void
    Simpler<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simpler model");
        {
          EquationOfState::LinearizedIncompressible<dim>::declare_parameters (prm);

          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          Rheology::ConstantViscosity::declare_parameters(prm,5e24);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Simpler<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simpler model");
        {
          equation_of_state.parse_parameters (prm);

          k_value                    = prm.get_double ("Thermal conductivity");

          constant_rheology.parse_parameters(prm);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::temperature;
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
    ASPECT_REGISTER_MATERIAL_MODEL(Simpler,
                                   "simpler",
                                   "A material model that has constant values "
                                   "except for density, which depends linearly on temperature: "
                                   "\\begin{align}"
                                   "  \\rho(p,T) &= \\left(1-\\alpha (T-T_0)\\right)\\rho_0."
                                   "\\end{align}"
                                   "\n\n"
                                   "\\note{This material model fills the role the ``simple'' material "
                                   "model was originally intended to fill, before the latter acquired "
                                   "all sorts of complicated temperature and compositional dependencies.}")
  }
}
