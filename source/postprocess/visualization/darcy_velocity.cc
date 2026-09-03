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


#include <aspect/melt.h>
#include <aspect/postprocess/visualization/darcy_velocity.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      DarcyVelocity<dim>::DarcyVelocity()
        : DataPostprocessorVector<dim>("darcy_velocity", update_values | update_quadrature_points | update_gradients)
        , Interface<dim>()
      {}



      template <int dim>
      std::string
      DarcyVelocity<dim>::get_physical_units() const
      {
        if (this->convert_output_to_years())
          return "m/year";
        else
          return "m/s";
      }



      template <int dim>
      void
      DarcyVelocity<dim>::evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                                std::vector<Vector<double>>                &computed_quantities) const
      {
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("The 'darcy velocity' postprocessor requires a compositional "
                               "field called porosity."));

        const unsigned int porosity_idx = this->introspection().find_composition_type(CompositionalFieldDescription::porosity);

        const double velocity_scaling_factor = this->convert_output_to_years() ? year_in_seconds : 1.0;

        MaterialModel::MaterialModelInputs<dim>  in(input_data, this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(in.n_evaluation_points(), this->n_compositional_fields());
        MeltHandler<dim>::create_material_model_outputs(out);
        this->get_material_model().evaluate(in, out);
        const std::shared_ptr<const MaterialModel::MeltOutputs<dim>> fluid_out =
          out.template get_additional_output_object<MaterialModel::MeltOutputs<dim>>();
        AssertThrow(std::isfinite(fluid_out->fluid_viscosities[0]),
                    ExcMessage("To compute the Darcy velocity the material model needs to provide the melt material model "
                               "outputs. At least the fluid viscosity was not computed, or is not a number."));

        for (unsigned int q = 0; q < in.n_evaluation_points(); ++q)
          {
            const Tensor<1, dim> gravity        = this->get_gravity_model().gravity_vector(in.position[q]);
            const Tensor<1, dim> solid_velocity = in.velocity[q];

            Tensor<1, dim> darcy_velocity = aspect::Utilities::calculate_approximate_darcy_velocity(
              in, out, fluid_out, solid_velocity, gravity, porosity_idx, q, this->get_parameters().use_pressure_gradient_for_darcy_field);
            for (unsigned int k = 0; k < dim; ++k)
              computed_quantities[q](k) = darcy_velocity[k] * velocity_scaling_factor;
          }
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(DarcyVelocity,
                                                  "darcy velocity",
                                                  "A visualization output object that outputs the Darcy velocity "
                                                  "vector. This postprocessor requires a compositional field named "
                                                  "'porosity'."
                                                  "\n\n"
                                                  "Physical units: $\\frac{\\text{m}}{\\text{s}}$ or "
                                                  "$\\frac{\\text{m}}{\\text{year}}$, depending on settings in the input file.")
    }
  }
}
