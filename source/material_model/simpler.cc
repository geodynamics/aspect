// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
    void
    Simpler<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // The Simpler model does not depend on composition
      EquationOfStateOutputs<dim> eos_outputs (1);

      thermal_conductivity.evaluate(in,out);

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          equation_of_state.evaluate(in, i, eos_outputs);

          out.viscosities[i] = constant_rheology.compute_viscosity();
          out.densities[i] = eos_outputs.densities[0];
          out.thermal_expansion_coefficients[i] = eos_outputs.thermal_expansion_coefficients[0];
          out.specific_heat[i] = eos_outputs.specific_heat_capacities[0];
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
          ThermalConductivity::Constant<dim>::declare_parameters (prm);
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
          thermal_conductivity.parse_parameters (prm);
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
                                   "$ \\rho(p,T) = \\left(1-\\alpha (T-T_0)\\right)\\rho_0.$ "
                                   "\n\n"
                                   "Note that this material model fills the role the ``simple'' material "
                                   "model was originally intended to fill, before the latter acquired "
                                   "all sorts of complicated temperature and compositional dependencies. ")
  }
}
