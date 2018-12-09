/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
      base_model_add->initialize();
    }


    template <int dim>
    void
    Additive<dim>::update()
    {
      base_model_add->update();
    }


    template <int dim>
    void
    Additive<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                            typename Interface<dim>::MaterialModelOutputs &out) const
    {
      const unsigned int density_idx = this->introspection().compositional_index_for_name("density_term");
      const unsigned int viscosity_idx = this->introspection().compositional_index_for_name("viscosity_factor");

      base_model_add -> evaluate(in,out);

      for (unsigned int i=0; i<in.position.size(); ++i)
        {
          out.densities[i] += in.composition[i][density_idx];
          out.viscosities[i] *= in.composition[i][viscosity_idx];

//          std::cout << "dens " << out.densities[i] << " visc " << out.viscosities[i] << std::endl;
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
          prm.declare_entry("Base model additive","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by an"
                            "averaging operation. Valid values for this parameter "
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
          Assert( prm.get("Base model additive") != "additive",
                  ExcMessage("You may not use ``additive'' as the base model for "
                             "a averaging model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model_add.reset(create_material_model<dim>(prm.get("Base model additive")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model_add.get()))
            sim->initialize_simulator (this->get_simulator());

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for averaging, it is essential to parse
      parameters related to the base model. */
      base_model_add->parse_parameters(prm);
      this->model_dependence = base_model_add->get_model_dependence();

    }


    template <int dim>
    bool
    Additive<dim>::
    is_compressible () const
    {
      return base_model_add->is_compressible();
    }

    template <int dim>
    double
    Additive<dim>::
    reference_viscosity() const
    {
      return base_model_add->reference_viscosity();
    }

//    template <int dim>
//    double
//    Additive<dim>::
//    reference_density() const
//    {
//      return base_model->reference_density();
//    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(Additive,
                                   "additive",
                                   " Explanation ")
  }
}

