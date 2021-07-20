/*
  Copyright (C) 2019 - 2021 by the authors of the ASPECT code.
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

#include <array>
#include <utility>
#include <limits>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/simulator.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/parsed_function.h>
#include <aspect/heating_model/interface.h>
#include <aspect/material_model/visco_plastic.h>

/* Head file for injection term*/
namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * This material model takes any other material model as a base model,
     * and adds additional material model outputs defining a dike injection
     * region of magma via a dilation term applied to the Stokes equations
     * that can be defined as a function depending on position and time in
     * the input file.
     *
     * The method is described in the following paper:
     * @code
     * @article{theissen2011coupled,
     *   title={Coupled mechanical and hydrothermal modeling of crustal accretion at intermediate to fast spreading ridges},
     *   author={Theissen-Krah, Sonja and Iyer, Karthik and R{\"u}pke, Lars H and Morgan, Jason Phipps},
     *   journal={Earth and Planetary Science Letters},
     *   volume={311},
     *   number={3-4},
     *   pages={275--286},
     *   year={2011},
     *   publisher={Elsevier}
     * }
     * @endcode
     *
     * @ingroup MaterialModels
     */

    template <int dim>
    class PrescribedDilation : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the model at the beginning of the run.
         */
        void initialize() override;

        /**
         * Update the base model and dilation function at the beginning
         * of each timestep.
         */
        void update() override;

        /**
         * Function to compute the material properties in @p out given
         * the inputs in @p in.
         */
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * Declare the parameters through input files.
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters through the input file
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Indicate whether material is compressible only based on
         * the base model.
         */
        bool is_compressible () const override;

        /**
         * Method to calculate reference viscosity.
         */
        double reference_viscosity () const override;

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:
        /**
         * Parsed function that specifies the region and amount of
         * material that is injected into the model.
         */
        Functions::ParsedFunction<dim> injection_function;

        /**
         * Pointer to the material model used as the base model.
         */
        std::shared_ptr<MaterialModel::Interface<dim> > base_model;
    };
  }
}

/* Head file for latent heat term*/
namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements the latent heat released during crystallization
     * of the melt lens and heating by melt injection into the model. It takes
     * the amount of material added on the right-hand side of the Stokes equations
     * and adds the corresponding heating term to the energy equation (considering
     * the latent heat of crystallization and the different temeprature of the
     * injected melt).
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class LatentHeatInjection : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Compute the heating model outputs for this class.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const override;
        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void parse_parameters (ParameterHandler &prm) override;
        /**
         * @}
         */

      private:
        /**
         * Properties of injected material.
         */
        double latent_heat_of_crystallization;
        double temperature_of_injected_melt;
    };
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    PrescribedDilation<dim>::initialize()
    {
      base_model->initialize();
    }

    template <int dim>
    void
    PrescribedDilation<dim>::update()
    {
      base_model->update();
    }

    template <int dim>
    void
    PrescribedDilation<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                      typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // fill variable out with the results form the base material model
      // The base model may have additional outputs such as frictional angle
      // in visco-plastic model, so we need to copy material properties first.
      typename Interface<dim>::MaterialModelOutputs base_output(out.n_evaluation_points(),
                                                                this->introspection().n_compositional_fields);

      // Move the additional outputs to base_output so that our models can fill
      // them if desired:
      base_output.move_additional_outputs_from(out);
      base_model->evaluate(in, base_output);

      // Copy required properties
      out.viscosities = base_output.viscosities;
      out.densities = base_output.densities;
      out.thermal_expansion_coefficients = base_output.thermal_expansion_coefficients;
      out.specific_heat = base_output.specific_heat;
      out.thermal_conductivities = base_output.thermal_conductivities;
      out.compressibilities = base_output.compressibilities;
      out.entropy_derivative_pressure = base_output.entropy_derivative_pressure;
      out.entropy_derivative_temperature = base_output.entropy_derivative_temperature;
      out.reaction_terms = base_output.reaction_terms;

      // Finally, we move the additional outputs back into place:
      out.move_additional_outputs_from(base_output);

      // Start to add the additional RHS terms to Stokes equations.
      MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force = (this->get_parameters().enable_additional_stokes_rhs)
               ? out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()
               : nullptr;

      MaterialModel::PrescribedPlasticDilation<dim>
      *prescribed_dilation = (this->get_parameters().enable_prescribed_dilation)
                             ? out.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim> >()
                             : nullptr;

      for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
        {
          // If "enable Enable additional Stokes RHS " is on, then the
          // additional rhs term is added into the rhs of the mass eq.
          if (force != nullptr)
            {
              for (unsigned int d=0; d < dim; ++d)
                force->rhs_u[i][d] = 0;

              // we get time passed as seconds (always) but may want
              // to reinterpret it in years.
              if (this->convert_output_to_years())
                force->rhs_p[i] = -1.0 * injection_function.value(in.position[i]) / year_in_seconds;
              else
                force->rhs_p[i] = -1.0 * injection_function.value(in.position[i]);
            }

          // If "Enable prescribed dilation" is on, then the injection rate
          // R (m/s or m/yr) is added to the rhs of the mass eq., i.e.,
          // -div(u,q) = -(R, q).
          // Meanwhile, the term - 2.0 / 3.0 * eta * (R, div v) is added to
          // the rhs of the momentum eq. (if the model is incompressible),
          // otherwise this term is already present on the left side.
          if (prescribed_dilation != nullptr)
            {
              if (this->convert_output_to_years())
                prescribed_dilation->dilation[i] = injection_function.value(in.position[i]) / year_in_seconds;
              else
                prescribed_dilation->dilation[i] = injection_function.value(in.position[i]);
            }

          // Mimic the oceanic crust generation/accretion in the dike.
          // We use the compositional reaction term to implement this process.
          const std::vector<double> &composition = in.composition[i];
          for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
            {
              // Lookup the injection area
              if (injection_function.value(in.position[i]) != 0.0)
                {
                  AssertThrow(this->introspection().compositional_name_exists("oceanic_crust"),
                              ExcMessage("Material model prescribed dilation only works if there "
                                         "is a compositional field called 'oceanic_crust'."));

                  if (c == this->introspection().compositional_index_for_name("oceanic_crust"))
                    out.reaction_terms[i][c] = -1.0 * composition[c] + 1.0;
                  else
                    out.reaction_terms[i][c] = -1.0 * composition[c];
                }
            }
        }

      // Currently, the dike injection process assumes that the nonlinear solver
      // scheme does a single advection iteration.
      // TODO: allow iterated advection.
      AssertThrow((this->get_parameters().nonlinear_solver ==
                   Parameters<dim>::NonlinearSolver::single_Advection_single_Stokes
                   ||
                   this->get_parameters().nonlinear_solver ==
                   Parameters<dim>::NonlinearSolver::single_Advection_iterated_Stokes
                   ||
                   this->get_parameters().nonlinear_solver ==
                   Parameters<dim>::NonlinearSolver::single_Advection_iterated_Newton_Stokes
                   ||
                   this->get_parameters().nonlinear_solver ==
                   Parameters<dim>::NonlinearSolver::single_Advection_iterated_defect_correction_Stokes),
                  ExcMessage("The material model will only work with the nonlinear "
                             "solver schemes 'single Advection, single Stokes', "
                             "'single Advection, iterated Stokes', "
                             "'single Advection, iterated Newton Stokes', and "
                             "'single Advection, iterated defect correction Stokes' "
                             "when prescribed dilation is enabled."));

    }

    template <int dim>
    void
    PrescribedDilation<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Prescribed dilation");
        {
          prm.declare_entry("Base model","simple",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by an "
                            "averaging operation. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.enter_subsection("Injection function");
          {
            Functions::ParsedFunction<dim>::declare_parameters(prm,1);
            prm.declare_entry("Function expression","0.0");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

    }

    template <int dim>
    void
    PrescribedDilation<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Prescribed dilation");
        {
          AssertThrow( prm.get("Base model") != "prescribed dilation",
                       ExcMessage("You may not use ''prescribed dilation'' as the base model for itself.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          prm.enter_subsection("Injection function");
          {
            try
              {
                injection_function.parse_parameters(prm);
              }
            catch (...)
              {
                std::cerr << "FunctionParser failed to parse\n"
                          << "\t Injection function\n"
                          << "with expression \n"
                          << "\t' " << prm.get("Function expression") << "'";
                throw;
              }
            prm.leave_subsection();
          }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // After parsing the parameters for averaging, it is essential
      // to parse parameters related to the base model.
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();

    }

    template <int dim>
    bool
    PrescribedDilation<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
    double
    PrescribedDilation<dim>::
    reference_viscosity() const
    {
      return true;
    }

    template <int dim>
    void
    PrescribedDilation<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // The base model may have additional outputs, so we need to copy
      // these additional outputs.
      base_model->create_additional_named_outputs(out);

      // Because we use the force outputs in the heating model, we always
      // have to attach them, not only in the places where the RHS of the
      // Stokes system is computed.
      const unsigned int n_points = out.n_evaluation_points();

      //Stokes additional RHS
      if (this->get_parameters().enable_additional_stokes_rhs
          && out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >() == nullptr)
        {
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>> (n_points));
        }

      AssertThrow(!this->get_parameters().enable_additional_stokes_rhs
                  ||
                  out.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()->rhs_u.size()
                  == n_points, ExcInternalError());

      // Prescribed dilation
      if (this->get_parameters().enable_prescribed_dilation
          && out.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim>>() == nullptr)
        {
          out.additional_outputs.push_back(
            std_cxx14::make_unique<MaterialModel::PrescribedPlasticDilation<dim>> (n_points));
        }

      AssertThrow(!this->get_parameters().enable_prescribed_dilation
                  ||
                  out.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim> >()->dilation.size()
                  == n_points, ExcInternalError());

    }
  }
}

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    LatentHeatInjection<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      AssertThrow(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
                  ExcMessage ("Heating outputs need to have the same number of entries as the material "
                              "model inputs."));

      const MaterialModel::AdditionalMaterialOutputsStokesRHS<dim>
      *force =
        (this->get_parameters().enable_additional_stokes_rhs)
        ? material_model_outputs.template get_additional_output<MaterialModel::AdditionalMaterialOutputsStokesRHS<dim> >()
        : nullptr;

      const MaterialModel::PrescribedPlasticDilation<dim>
      *prescribed_dilation =
        (this->get_parameters().enable_prescribed_dilation)
        ? material_model_outputs.template get_additional_output<MaterialModel::PrescribedPlasticDilation<dim> >()
        : nullptr;

      // Add the latent heat source term corresponding to prescribed dilation
      // terms in Stokes equations to the rhs of energy conservation equation.
      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = 0.0;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
          heating_model_outputs.rates_of_temperature_change[q] = 0.0;

          if (force != nullptr)
            heating_model_outputs.heating_source_terms[q] = -1.0 * force->rhs_p[q] * (latent_heat_of_crystallization + (temperature_of_injected_melt - material_model_inputs.temperature[q]) * material_model_outputs.densities[q] * material_model_outputs.specific_heat[q]);

          if (prescribed_dilation != nullptr)
            heating_model_outputs.heating_source_terms[q] = prescribed_dilation->dilation[q] * (latent_heat_of_crystallization + (temperature_of_injected_melt - material_model_inputs.temperature[q]) * material_model_outputs.densities[q] * material_model_outputs.specific_heat[q]);

        }
    }

    template <int dim>
    void
    LatentHeatInjection<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat injection");
        {
          prm.declare_entry ("Latent heat of crystallization", "1.1e9",
                             Patterns::Double(0),
                             "The latent heat of crystallization that is released when material "
                             "is injected into the model. "
                             "Units: J/m$^3$.");
          prm.declare_entry ("Temperature of injected melt", "1573",
                             Patterns::Double(0),
                             "The temperature of the material injected into the model. "
                             "Units: K.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    LatentHeatInjection<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const
    {
      this->get_material_model().create_additional_named_outputs(material_model_outputs);
    }

    template <int dim>
    void
    LatentHeatInjection<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Latent heat injection");
        {
          latent_heat_of_crystallization = prm.get_double ("Latent heat of crystallization");
          temperature_of_injected_melt = prm.get_double ("Temperature of injected melt");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PrescribedDilation,
                                   "prescribed dilation",
                                   "The material model uses a ``Base model'' from which material properties are "
                                   "derived. It then adds source terms in the Stokes equations "
                                   "that describe a dike injection of melt to the model. ")
  }

  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(LatentHeatInjection,
                                  "latent heat injection",
                                  "Latent heat releases due to the injection of melt into the model. "
                                  "This heating model takes the source term added to the Stokes "
                                  "equation and adds the corresponding source term to the energy "
                                  "equation. This source term includes both the effect of latent "
                                  "heat release upon crystallization and the fact that inected "
                                  "material might have a diffeent temperature.")
  }
}
