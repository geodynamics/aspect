#include <deal.II/base/function_lib.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <deal.II/base/parameter_handler.h>

#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class ExponentialDecay : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
        static void declare_parameters (ParameterHandler &prm);
        virtual void parse_parameters(ParameterHandler &prm);

        virtual bool is_compressible () const;
        virtual double reference_viscosity () const;

      private:

        /**
         * Pointer to the material model used as the base model
         */
        std_cxx11::shared_ptr<MaterialModel::Interface<dim> > base_model;

        /**
         * Parameters determining the initial decay rate.
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
    double
    ExponentialDecay<dim>::
    reference_viscosity() const
    {
      return base_model->reference_viscosity();
    }


    template <int dim>
    void
    ExponentialDecay<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      AssertThrow(this->introspection().n_compositional_fields == 1,
                  ExcMessage("Exponential decay model needs exactly one compositiona field."));

      // Fill material model outputs using the base model.
      base_model->evaluate(in,out);
      const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

      for (unsigned int q=0; q < in.position.size(); ++q)
        for (unsigned int c=0; c < this->introspection().n_compositional_fields; ++c)
          out.reaction_terms[q][c] = 0.0;

      // fill melt reaction rates if they exist
      ReactionRateOutputs<dim> *reaction_out = out.template get_additional_output<ReactionRateOutputs<dim> >();

      if (reaction_out != NULL)
        {
          for (unsigned int q=0; q < in.position.size(); ++q)
            {
              // dx/dt = alpha * x + beta * x * y
              const double decay_constant = half_life > 0.0 ? log(2.0) / half_life : 0.0;
              reaction_out->reaction_rates[0][q] = - decay_constant / time_scale * in.composition[q][0];
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
          base_model.reset(create_material_model<dim>(prm.get("Base model")));
          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
            sim->initialize_simulator (this->get_simulator());

          half_life              = prm.get_double ("Half life");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      /* After parsing the parameters for depth dependent, it is essential to parse
      parameters related to the base model. */
      base_model->parse_parameters(prm);
      this->model_dependence = base_model->get_model_dependence();
    }
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
                                   "function that can be chosen as an input parameter "
                                   "(units: 1/s or 1/yr).")

  }
}

