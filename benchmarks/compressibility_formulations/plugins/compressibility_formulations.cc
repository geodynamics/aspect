/*
  Copyright (C) 2017 - 2024 by the authors of the ASPECT code.

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


#include <aspect/simulator/assemblers/interface.h>
#include <aspect/simulator/assemblers/stokes.h>
#include <aspect/material_model/simple_compressible.h>
#include <aspect/simulator_access.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    /**
     * A material model that is identical to the simple compressible model,
     * except that the density is tracked in a compositional field of type
     * 'density' using the prescribed field advection method. It also
     * allows some modification to the density and thermal expansivity
     * calculation for the compressibility benchmarks.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class CompressibilityFormulations : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void initialize() override;

        void update() override;

        bool is_compressible () const override;

        void
        evaluate (const MaterialModelInputs<dim> &in,
                  MaterialModelOutputs<dim> &out) const override;

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

        static void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Pointer to the material model used as the base model
         */
        std::shared_ptr<MaterialModel::Interface<dim>> base_model;

        /**
         * The reference density
         */
        double reference_rho;

        /**
         * The constant thermal expansivity
         */
        double thermal_alpha;
        bool use_exponential_alpha;

        /**
         * The constant compressibility.
         */
        double reference_compressibility;

        /**
         * Use the adiabatic instead of the full pressure for the
         * density calculation.
         */
        bool use_adiabatic_pressure_for_density;
    };



    template <int dim>
    void
    CompressibilityFormulations<dim>::initialize()
    {
      base_model->initialize();
    }



    template <int dim>
    void
    CompressibilityFormulations<dim>::update()
    {
      base_model->update();
    }



    template <int dim>
    bool
    CompressibilityFormulations<dim>::
    is_compressible () const
    {
      return this->get_parameters().formulation_mass_conservation != Parameters<dim>::Formulation::MassConservation::incompressible;
    }



    template <int dim>
    void
    CompressibilityFormulations<dim>::evaluate(const MaterialModelInputs<dim> &in,
                                               MaterialModelOutputs<dim> &out) const
    {
      base_model->evaluate(in,out);

      const unsigned int density_field_index = this->introspection().find_composition_type(Parameters<dim>::CompositionalFieldDescription::density);

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          double pressure_for_density = in.pressure[i];

          // If using the PDA avoid using the real pressure for the density calculation to prevent
          // pressure waves from forming
          if (use_adiabatic_pressure_for_density)
            pressure_for_density = this->get_adiabatic_conditions().pressure(in.position[i]);

          double thermal_expansivity = thermal_alpha;
          if (use_exponential_alpha)
            {
              thermal_expansivity *= std::exp(-1.117979323e-11*pressure_for_density);
              out.thermal_expansion_coefficients[i] = thermal_expansivity;
            }

          // Use a thermodynamically consistent thermal expansivity
          // (if alpha is constant, rho needs to depend exponentially on it, not linearly)
          out.densities[i] = reference_rho * std::exp(reference_compressibility * (pressure_for_density - this->get_surface_pressure()) -
                                                      thermal_expansivity * (in.temperature[i] - this->get_adiabatic_surface_temperature()));

          // For the ICA, we have to provide the isentropic rather than isothermal compressibility
          // for the mass conservation equation.
          if (this->get_parameters().formulation_mass_conservation ==
              Parameters<dim>::Formulation::MassConservation::isentropic_compression)
            {
              out.compressibilities[i] = reference_compressibility - Utilities::fixed_power<2>(thermal_expansivity) * in.temperature[i]
                                         / (out.densities[i] * out.specific_heat[i]);
            }
        }

      // prescribe the density field to the current value of the density. The actual projection
      // only happens inside Simulator<dim>::interpolate_material_output_into_compositional_field,
      // this just sets the correct term the field will be set to.
      if (PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>())
        for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
          prescribed_field_out->prescribed_field_outputs[i][density_field_index] = out.densities[i];
    }



    template <int dim>
    void
    CompressibilityFormulations<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PrescribedFieldOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::unique_ptr<MaterialModel::AdditionalMaterialOutputs<dim>>
            (new MaterialModel::PrescribedFieldOutputs<dim> (n_points, this->n_compositional_fields())));
        }
    }



    template <int dim>
    void
    CompressibilityFormulations<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Projected density model");
        {
          prm.declare_entry ("Use adiabatic pressure for density", "false",
                             Patterns::Bool(),
                             "");
          prm.declare_entry ("Use exponential thermal expansivity", "false",
                             Patterns::Bool(),
                             "");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      SimpleCompressible<dim>::declare_parameters (prm);
    }



    template <int dim>
    void
    CompressibilityFormulations<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Projected density model");
        {
          use_adiabatic_pressure_for_density = prm.get_bool("Use adiabatic pressure for density");
          use_exponential_alpha = prm.get_bool("Use exponential thermal expansivity");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();


      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Simple compressible model");
        {
          reference_rho              = prm.get_double ("Reference density");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          reference_compressibility  = prm.get_double ("Reference compressibility");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // create the base model and initialize its SimulatorAccess base
      // class; it will get a chance to read its parameters below after we
      // leave the current section
      base_model = create_material_model<dim>("simple compressible");
      if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
        sim->initialize_simulator (this->get_simulator());
      base_model->parse_parameters(prm);
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(CompressibilityFormulations,
                                   "compressibility formulations",
                                   "A material model that is identical to the simple compressible model, "
                                   "except that the density is tracked in a compositional field of type "
                                   "'density' using the prescribed field advection method. It also "
                                   "allows some modification to the density and thermal expansivity "
                                   "calculation for the compressibility benchmarks.")
  }
}
