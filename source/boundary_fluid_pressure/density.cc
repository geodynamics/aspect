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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/boundary_fluid_pressure/density.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/melt.h>
#include <utility>
#include <limits>


namespace aspect
{
  namespace BoundaryFluidPressure
  {

    template <int dim>
    void
    Density<dim>::
    fluid_pressure_gradient (
      const types::boundary_id /*boundary_indicator*/,
      const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
      const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
      const std::vector<Tensor<1,dim> > &normal_vectors,
      std::vector<double> &fluid_pressure_gradient_outputs
    ) const
    {
      const MaterialModel::MeltOutputs<dim> *melt_outputs = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim> >();
      Assert(melt_outputs!=NULL, ExcMessage("Error, MeltOutputs are missing in fluid_pressure_gradient()"));
      for (unsigned int q=0; q<fluid_pressure_gradient_outputs.size(); ++q)
        {
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(material_model_inputs.position[q]);

          switch (density_formulation)
            {
              case DensityFormulation::solid_density:
              {
                fluid_pressure_gradient_outputs[q] = (material_model_outputs.densities[q] * gravity) * normal_vectors[q];
                break;
              }

              case DensityFormulation::fluid_density:
              {
                fluid_pressure_gradient_outputs[q] = (melt_outputs->fluid_densities[q] * gravity) * normal_vectors[q];
                break;
              }

              default:
                Assert (false, ExcNotImplemented());
            }
        }
    }

    template <int dim>
    void
    Density<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary fluid pressure model");
      {
        prm.enter_subsection("Density");
        {
          prm.declare_entry ("Density formulation", "solid density",
                             Patterns::Selection ("solid density|fluid density"),
                             "The density formulation used to compute the fluid pressure gradient "
                             "at the model boundary."
                             "\n\n"
                             "`solid density' prescribes the gradient of the fluid pressure as "
                             "solid density times gravity (which is the lithostatic "
                             "pressure) and leads to approximately the same pressure in "
                             "the melt as in the solid, so that fluid is only flowing "
                             "in or out due to differences in dynamic pressure."
                             "\n\n"
                             "`fluid density' prescribes the gradient of the fluid pressure as "
                             "fluid density times gravity and causes melt to flow in "
                             "with the same velocity as inflowing solid material, "
                             "or no melt flowing in or out if the solid velocity "
                             "normal to the boundary is zero.");
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Density<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Boundary fluid pressure model");
      {
        prm.enter_subsection("Density");
        {
          if (prm.get ("Density formulation") == "solid density")
            density_formulation = DensityFormulation::solid_density;
          else if (prm.get ("Density formulation") == "fluid density")
            density_formulation = DensityFormulation::fluid_density;
          else
            AssertThrow (false, ExcNotImplemented());
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace BoundaryFluidPressure
  {
    ASPECT_REGISTER_BOUNDARY_FLUID_PRESSURE_MODEL(Density,
                                                  "density",
                                                  "A plugin that prescribes the fluid pressure gradient at "
                                                  "the boundary based on fluid/solid density from the material "
                                                  "model.")
  }
}
