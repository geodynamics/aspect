/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#include <aspect/fluid_pressure_boundary_conditions/density.h>
#include <aspect/melt.h>
#include <utility>
#include <limits>


namespace aspect
{
  namespace FluidPressureBoundaryConditions
  {

    template <int dim>
    void
    Density<dim>::
    fluid_pressure_gradient (
      const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
      const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
      std::vector<Tensor<1,dim> > &output
    ) const
    {
      const MaterialModel::MeltOutputs<dim> *melt_outputs = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
      Assert(melt_outputs!=NULL, ExcMessage("Error, MeltOutputs are missing in fluid_pressure_gradient()"));
      for (unsigned int q=0; q<output.size(); ++q)
        {
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(material_model_inputs.position[q]);
          output[q] = (((include_rho_s) ? material_model_outputs.densities[q] : 0.0)
                       + ((include_rho_f) ? melt_outputs->fluid_densities[q] : 0.0))
                      * gravity;
        }
    }

    template <int dim>
    void
    Density<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Fluid Pressure Boundary Condition");
      {
        prm.enter_subsection("Density");
        {
          prm.declare_entry ("Include Fluid Density", "false", Patterns::Bool(),
                             "Prescribing the gradient of the fluid pressure as "
                             "fluid density times gravity causes melt to flow in "
                             "with the same velocity as inflowing solid material, "
                             "or no melt flowing in or out if the solid velocity "
                             "normal to the boundary is zero.");
          prm.declare_entry ("Include Solid Density", "true", Patterns::Bool(),
                             "Prescribing the gradient of the fluid pressure as "
                             "solid density times gravity (which is the lithostatic "
                             "pressure) leads to approximately the same pressure in "
                             "the melt as in the solid, so that fluid is only flowing "
                             "in or out due to differences in dynamic pressure. ");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Density<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Fluid Pressure Boundary Condition");
      {
        prm.enter_subsection("Density");
        {
          include_rho_f = prm.get_bool("Include Fluid Density");
          include_rho_s = prm.get_bool("Include Solid Density");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    Density<dim>::initialize()
    {
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace FluidPressureBoundaryConditions
  {
    ASPECT_REGISTER_FLUID_PRESSURE_BOUNDARY_CONDITIONS(Density,
                                                       "density",
                                                       "A plugin that prescribes the fluid pressure gradient at "
                                                       "the boundary based on fluid/solid density from the material "
                                                       "model. ")
  }
}
