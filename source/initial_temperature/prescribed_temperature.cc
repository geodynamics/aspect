/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/initial_temperature/prescribed_temperature.h>

#include <aspect/initial_composition/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>


namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    void
    PrescribedTemperature<dim>::initialize()
    {
      // Make sure we keep track of the initial composition manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_composition = this->get_initial_composition_manager_pointer();
    }


    template <int dim>
    double
    PrescribedTemperature<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // Evaluate the material model to get the temperature
      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());
      in.requested_properties = MaterialModel::MaterialProperties::additional_outputs;

      in.position[0] = position;
      in.temperature[0] = this->get_adiabatic_conditions().temperature(position);
      in.pressure[0] = this->get_adiabatic_conditions().pressure(position);
      in.velocity[0] = Tensor<1,dim> ();

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        in.composition[0][c] = initial_composition->initial_composition(position, c);

      in.strain_rate.resize(0);

      this->get_material_model().create_additional_named_outputs(out);
      this->get_material_model().evaluate(in, out);

      // set up variable to interpolate prescribed field outputs onto temperature field
      double temperature = in.temperature[0];
      if (MaterialModel::PrescribedTemperatureOutputs<dim> *prescribed_temperature_out
          = out.template get_additional_output<MaterialModel::PrescribedTemperatureOutputs<dim>>())
        {
          temperature = prescribed_temperature_out->prescribed_temperature_outputs[0];
        }
      else
        AssertThrow (false,
                     ExcMessage ("The material model needs to compute precribed temperature outputs, "
                                 "otherwise the initial temperature can not be set to the prescribed "
                                 "temperature."));

      return temperature;
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(PrescribedTemperature,
                                              "prescribed temperature",
                                              "This model fixes the initial temperature to the prescribed "
                                              "temperature outputs computed by the material model. This only "
                                              "works if the material model implements prescribed temperature "
                                              "outputs.")
  }
}
