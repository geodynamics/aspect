/*
  Copyright (C) 2021 by the authors of the ASPECT code.

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

#include <aspect/material_model/additive.h>

using namespace dealii;

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    Additive<dim>::initialize()
    {

      // This material model is only useful when using and adjoint Stokes solver scheme so throw
      // an exception if a different solver is chosen.
      Assert(this->get_parameters().nonlinear_solver == Parameters<dim>::NonlinearSolver::no_Advection_adjoint_Stokes,
             ExcMessage("Error: Material model Additive is only sensible with the no Advection, adjoint Stokes solver scheme."))

      // If the no Advection, adjoint Stokes solver scheme is chosen the parameter file is already checked such
      // that it has at least 2 compositional fields and that the compositional fields have the correct names.

      base_model->initialize();
    }



    template <int dim>
    void
    Additive<dim>::update()
    {
      base_model->update();
    }



    template <int dim>
    void
    Additive<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                            typename Interface<dim>::MaterialModelOutputs &out) const
    {

      const unsigned int density_idx = this->introspection().compositional_index_for_name("density_increment");
      const unsigned int viscosity_idx = this->introspection().compositional_index_for_name("viscosity_increment");

      base_model -> evaluate(in,out);

      for (unsigned int i=0; i<in.position.size(); ++i)
        {
          // add the values from the compositional field to the density / viscosity to update
          // the parameters. This allows material properties to change over several adjoint iterations.
          // Note that the viscosity is already multiplied with the background viscosity.
          out.densities[i] += in.composition[i][density_idx];
          out.viscosities[i] += in.composition[i][viscosity_idx];

        }
    } 



    template <int dim>
    void
    Additive<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Additive model");
        {
          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by an"
                            "addition operation. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Additive<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Additive model");
        {
          Assert( prm.get("Base model") != "additive",
                  ExcMessage("You may not use ``additive'' as the base model for "
                             "an additive model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // After parsing the parameters for averaging, it is essential to parse
      // parameters related to the base model.
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }

    template <int dim>
    bool
    Additive<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
    double
    Additive<dim>::
    reference_viscosity() const
    {
      return base_model->reference_viscosity();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Additive,
                                   "additive",
                                   "This material model is only sensible with the adjoint "
                                   "Stokes nonlinear solver scheme. "
                                   "The material model allows the density and viscosity "
                                   "to change over several iterations within the same "
                                   "timestep. The change in parameters is given by the "
                                   "the two compositional fields. This material model "
                                   "therefore needs to be run with at least 2 compositional fields. "
                                   "with the names 'density_increment' and 'viscosity_increment'. "
                                   "The compositional fields are updated to minimize "
                                   "the objective functional in the adjoint equations.")
  }
}
