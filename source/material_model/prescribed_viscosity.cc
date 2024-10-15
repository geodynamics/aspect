/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/material_model/prescribed_viscosity.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    PrescribedViscosity<dim>::initialize()
    {
      base_model->initialize();
    }



    template <int dim>
    void
    PrescribedViscosity<dim>::update()
    {
      base_model->update();

      // we get time passed as seconds (always) but may want
      // to reinterpret it in years
      const double time = this->convert_output_to_years() ? this->get_time() / year_in_seconds : this->get_time();
      prescribed_viscosity_indicator_function.set_time (time);
      prescribed_viscosity_function.set_time (time);
    }



    template <int dim>
    void
    PrescribedViscosity<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                       typename Interface<dim>::MaterialModelOutputs &out) const
    {
      base_model->evaluate(in,out);
      for (unsigned int i=0; i<out.n_evaluation_points(); ++i)
        {
          if (in.requests_property(MaterialProperties::viscosity))
            {
              const double indicator = prescribed_viscosity_indicator_function.value(in.position[i]);

              if (indicator > 0.5)
                {
                  out.viscosities[i] = prescribed_viscosity_function.value(in.position[i]);
                }
            }
        }
    }



    template <int dim>
    bool
    PrescribedViscosity<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }



    template <int dim>
    void
    PrescribedViscosity<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Prescribed viscosity");
        {

          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by the prescribed "
                            "viscosity material model. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");

          prm.enter_subsection ("Indicator function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 3);
          }
          prm.leave_subsection ();

          prm.enter_subsection ("Viscosity function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 3);
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    PrescribedViscosity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection ("Prescribed viscosity");
        {
          AssertThrow( prm.get("Base model") != "prescribed viscosity",
                       ExcMessage("You may not use ``prescribed viscosity'' as the base model for "
                                  "a prescribed viscosity model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model = create_material_model<dim>(prm.get("Base model"));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          prm.enter_subsection("Indicator function");
          {
            try
              {
                prescribed_viscosity_indicator_function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed viscosity.Indicator function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'";
                throw;
              }
          }
          prm.leave_subsection();

          prm.enter_subsection("Viscosity function");
          {
            try
              {
                prescribed_viscosity_function.parse_parameters (prm);
              }
            catch (...)
              {
                std::cerr << "ERROR: FunctionParser failed to parse\n"
                          << "\t'Prescribed viscosity.Viscosity function'\n"
                          << "with expression\n"
                          << "\t'" << prm.get("Function expression") << "'";
                throw;
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();

      /* After parsing the parameters for prescribed viscosity, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }



    template <int dim>
    void
    PrescribedViscosity<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      base_model->create_additional_named_outputs(out);
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PrescribedViscosity,
                                   "prescribed viscosity",
                                   "A material model that applies a viscosity to a ''base model'' chosen from any of "
                                   "the other available material models. This prescribed viscosity material model "
                                   "allows the user to specify a function which describes where the viscosity should be "
                                   "prescribed and a second function which describes the viscosity in that region. "
                                   "This material model requires a base model which prescribes the viscosity and the "
                                   "other material parameters in the rest of the model."
                                  )
  }
}
