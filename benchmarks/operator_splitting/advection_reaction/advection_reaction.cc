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
#include <deal.II/base/function_lib.h>
#include <aspect/material_model/interface.h>
#include <aspect/boundary_composition/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/vector_tools.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    class ExponentialDecay : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;
        static void declare_parameters (ParameterHandler &prm);
        void parse_parameters(ParameterHandler &prm) override;

        bool is_compressible () const override;

        void create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        /**
         * Pointer to the material model used as the base model
         */
        std::unique_ptr<MaterialModel::Interface<dim>> base_model;

        /**
         * Parameters determining the initial decay rate.
         */
        double half_life;
    };
  }

  namespace HeatingModel
  {
    template <int dim>
    class ExponentialDecayHeating : public HeatingModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const override;

        static void declare_parameters (ParameterHandler &prm);

        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Parameter determining the decay rate.
         */
        double half_life;
    };
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    bool
    ExponentialDecay<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }


    template <int dim>
    void
    ExponentialDecay<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      AssertThrow(this->introspection().n_compositional_fields == 1,
                  ExcMessage("Exponential decay model needs exactly one compositional field."));

      // Fill material model outputs using the base model.
      base_model->evaluate(in,out);
      const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

      for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
        for (unsigned int c=0; c < this->introspection().n_compositional_fields; ++c)
          out.reaction_terms[q][c] = 0.0;

      // fill melt reaction rates if they exist
      ReactionRateOutputs<dim> *reaction_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

      if (reaction_out != nullptr)
        {
          for (unsigned int q=0; q < in.n_evaluation_points(); ++q)
            {
              // dC/dt = - z * lambda * C
              const double decay_constant = half_life > 0.0 ? std::log(2.0) / half_life : 0.0;
              const double z = in.position[q](1);
              reaction_out->reaction_rates[q][0] = - decay_constant * z / time_scale * in.composition[q][0];
            }
        }
    }


    template <int dim>
    void
    ExponentialDecay<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Exponential decay");
        {
          prm.declare_entry("Base model","composition reaction",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a position "
                            "dependent melting rate. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");
          prm.declare_entry ("Half life", "0",
                             Patterns::Double (0),
                             "Time required for a compositional field to reduce to half its "
                             "initial value. Units: Years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ExponentialDecay<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Exponential decay");
        {
          AssertThrow( prm.get("Base model") != "reaction rate function",
                       ExcMessage("You may not use ``reaction rate function'' as the base model for "
                                  "a reaction rate function model.") );

          // create the base model and initialize its SimulatorAccess base
          // class; it will get a chance to read its parameters below after we
          // leave the current section
          base_model = create_material_model<dim>(prm.get("Base model"));
          if (auto s = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            s->initialize_simulator (this->get_simulator());

          half_life              = prm.get_double ("Half life");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // After parsing the parameters for the exponential decay material model,
      // also parse the parameters related to the base model.
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }


    template <int dim>
    void
    ExponentialDecay<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }


  namespace HeatingModel
  {
    template <int dim>
    void
    ExponentialDecayHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
              const MaterialModel::MaterialModelOutputs<dim> & /*out*/,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == in.n_evaluation_points(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

      for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
        {
          // dC/dt = - z * lambda * C
          const double decay_constant = half_life > 0.0 ? std::log(2.0) / half_life : 0.0;
          const double z = in.position[q](1);
          heating_model_outputs.rates_of_temperature_change[q] = - decay_constant * z / time_scale * in.composition[q][0];

          heating_model_outputs.heating_source_terms[q] = 0.0;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
        }
    }


    template <int dim>
    void
    ExponentialDecayHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Exponential decay heating");
        {
          prm.declare_entry ("Half life", "0",
                             Patterns::Double (0),
                             "Time required for the temperature to reduce to half of its "
                             "initial value. Units: Years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    ExponentialDecayHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Exponential decay heating");
        {
          half_life              = prm.get_double ("Half life");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }


  template <int dim>
  class RefFunction : public Function<dim>
  {
    public:
      RefFunction () : Function<dim>(dim+3) {}

      void vector_value (const Point<dim>   &p,
                         Vector<double>   &values) const override
      {
        double x = p(0);
        double z = p(1);
        double t = this->get_time();

        values[0] = 0.0; // velocity x
        values[1] = 0.0; // velocity z
        values[2] = 0.0; // pressure
        values[3] = std::sin(2*numbers::PI*(x-t*0.01))*std::exp(-std::log(2.0)/10.0*z*t); // temperature
        values[4] = std::sin(2*numbers::PI*(x-t*0.01))*std::exp(-std::log(2.0)/10.0*z*t); // composition
      }
  };



  template <int dim>
  class ExponentialDecayBoundary : public BoundaryComposition::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Return the boundary composition as a function of position.
       */
      double
      boundary_composition (const types::boundary_id boundary_indicator,
                            const Point<dim> &position,
                            const unsigned int compositional_field) const override;
  };


  template <int dim>
  double
  ExponentialDecayBoundary<dim>::
  boundary_composition (const types::boundary_id /*boundary_indicator*/,
                        const Point<dim> &p,
                        const unsigned int /*compositional_field*/) const
  {
    double x = p(0);
    double z = p(1);
    const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

    double t = this->get_time() / time_scale;

    return std::sin(2*numbers::PI*(x-t*0.01))*std::exp(-std::log(2.0)/10.0*z*t);
  }



  template <int dim>
  class ExponentialDecayBoundaryTemperature : public BoundaryTemperature::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      /**
       * Return the boundary composition as a function of position.
       */
      double
      boundary_temperature (const types::boundary_id boundary_indicator,
                            const Point<dim> &position) const override;

      double minimal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const override;

      double maximal_temperature (const std::set<types::boundary_id> &fixed_boundary_ids) const override;
  };


  template <int dim>
  double
  ExponentialDecayBoundaryTemperature<dim>::
  boundary_temperature (const types::boundary_id /*boundary_indicator*/,
                        const Point<dim> &p) const
  {
    double x = p(0);
    double z = p(1);
    const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

    double t = this->get_time() / time_scale;

    return std::sin(2*numbers::PI*(x-t*0.01))*std::exp(-std::log(2.0)/10.0*z*t);
  }

  template <int dim>
  double
  ExponentialDecayBoundaryTemperature<dim>::
  minimal_temperature (const std::set<types::boundary_id> &/*fixed_boundary_ids*/) const
  {
    return 0.0;
  }



  template <int dim>
  double
  ExponentialDecayBoundaryTemperature<dim>::
  maximal_temperature (const std::set<types::boundary_id> &/*fixed_boundary_ids*/) const
  {
    return 1.0;
  }



  /**
     * A postprocessor that evaluates the accuracy of the solution
     * by using the L2 norm.
     */
  template <int dim>
  class ExponentialDecayPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
  {
    public:
      ExponentialDecayPostprocessor();

      /**
       * Generate graphical output from the current solution.
       */
      std::pair<std::string,std::string>
      execute (TableHandler &statistics) override;

      double max_error;
      double max_error_T;
  };

  template <int dim>
  ExponentialDecayPostprocessor<dim>::ExponentialDecayPostprocessor ()
  {
    max_error = 0.0;
    max_error_T = 0.0;
  }

  template <int dim>
  std::pair<std::string,std::string>
  ExponentialDecayPostprocessor<dim>::execute (TableHandler & /*statistics*/)
  {
    AssertThrow (this->n_compositional_fields() == 1, ExcInternalError());

    RefFunction<dim> ref_func;
    ref_func.set_time(this->get_time());

    const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

    const unsigned int n_total_comp = this->introspection().n_components;

    Vector<float> cellwise_errors_composition (this->get_triangulation().n_active_cells());
    Vector<float> cellwise_errors_temperature (this->get_triangulation().n_active_cells());

    ComponentSelectFunction<dim> comp_C(this->introspection().component_indices.compositional_fields[0], n_total_comp);
    ComponentSelectFunction<dim> comp_T(this->introspection().component_indices.temperature, n_total_comp);

    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_composition,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_C);
    VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                       this->get_solution(),
                                       ref_func,
                                       cellwise_errors_temperature,
                                       quadrature_formula,
                                       VectorTools::L2_norm,
                                       &comp_T);

    const double current_error = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_composition, VectorTools::L2_norm);
    max_error = std::max(max_error, current_error);

    const double current_error_T = VectorTools::compute_global_error(this->get_triangulation(), cellwise_errors_temperature, VectorTools::L2_norm);
    max_error_T = std::max(max_error_T, current_error_T);

    std::ostringstream os;
    os << std::scientific
       << "time=" << this->get_time()
       << " ndofs= " << this->get_solution().size()
       << " C_L2_current= " << current_error
       << " C_L2_max= " << max_error
       << " T_L2_current= " << current_error_T
       << " T_L2_max= " << max_error_T
       ;

    return std::make_pair("Errors", os.str());
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ExponentialDecay,
                                   "exponential decay",
                                   "A material model that can be derived from any of the other "
                                   "material model and that will replace the reaction rate by a "
                                   "function that models exponential decay. The half life can be "
                                   "chosen as an input parameter.")
  }
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ExponentialDecayHeating,
                                  "exponential decay heating",
                                  "A heating model that will use a model for exponential decay as "
                                  "the heating reaction rate. The half life can be chosen as an "
                                  "input parameter.")

  }
  ASPECT_REGISTER_BOUNDARY_COMPOSITION_MODEL(ExponentialDecayBoundary,
                                             "ExponentialDecayBoundary",
                                             "A boundary composition model that prescribes the analytical "
                                             "solution for an initially sinusoidal field that is advected "
                                             "and decays over time.")

  ASPECT_REGISTER_BOUNDARY_TEMPERATURE_MODEL(ExponentialDecayBoundaryTemperature,
                                             "ExponentialDecayBoundary",
                                             "A boundary composition model that prescribes the analytical "
                                             "solution for an initially sinusoidal field that is advected "
                                             "and decays over time.")

  ASPECT_REGISTER_POSTPROCESSOR(ExponentialDecayPostprocessor,
                                "ExponentialDecayPostprocessor",
                                "A postprocessor that compares the solution "
                                "to the analytical solution for exponential decay.")
}
