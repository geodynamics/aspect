#include <deal.II/base/function_lib.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class MeltingRateFunction : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        MeltingRateFunction ();

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
         * A function object representing the melting rate.
         */
        Functions::ParsedFunction<dim> function;
    };

  }
}

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    MeltingRateFunction<dim>::MeltingRateFunction ()
    :
    function (1)
  {}


    template <int dim>
    bool
    MeltingRateFunction<dim>::
    is_compressible () const
    {
      return base_model->is_compressible();
    }

    template <int dim>
    double
    MeltingRateFunction<dim>::
    reference_viscosity() const
    {
      return base_model->reference_viscosity();
    }


    template <int dim>
    void
    MeltingRateFunction<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      // Fill material model outputs using the base model.
      base_model->evaluate(in,out);
      const double time_scale = this->convert_output_to_years() ? year_in_seconds : 1.0;

      // overwrite the melting rate
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
      const unsigned int peridotite_idx = this->introspection().compositional_index_for_name("peridotite");
      if (in.cell && this->get_timestep_number() > 0)
        {
          for (unsigned int q=0; q < in.position.size(); ++q)
            {
              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  if (c == porosity_idx)
                    out.reaction_terms[q][c] = function.value(in.position[q]) / time_scale * out.densities[q];
                  else if (c == peridotite_idx)
                    out.reaction_terms[q][c] = 0.0;
                }
              out.viscosities[q] *= (1.0 - in.composition[q][porosity_idx]);
          }
        }

      // fill melt outputs if they exist
//      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim> >();
//
//      if (melt_out != NULL)
//        {
//          for (unsigned int q=0; q < in.position.size(); ++q)
//            {
//              melt_out->compaction_viscosities[q] *= (1.0 - in.composition[q][porosity_idx]);
//            }
//        }
    }


    template <int dim>
    void
    MeltingRateFunction<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melting rate function");
        {
          prm.declare_entry("Base model","melt visco plastic",
                            Patterns::Selection(MaterialModel::get_valid_model_names_pattern<dim>()),
                            "The name of a material model that will be modified by a position "
                            "dependent melting rate. Valid values for this parameter "
                            "are the names of models that are also valid for the "
                            "``Material models/Model name'' parameter. See the documentation for "
                            "that for more information.");

          prm.enter_subsection("Function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltingRateFunction<dim>::
    parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
          prm.enter_subsection("Melting rate function");
          {
              AssertThrow( prm.get("Base model") != "melting rate function",
                           ExcMessage("You may not use ``depth dependent'' as the base model for "
                                      "a depth-dependent model.") );

              // create the base model and initialize its SimulatorAccess base
              // class; it will get a chance to read its parameters below after we
              // leave the current section
              base_model.reset(create_material_model<dim>(prm.get("Base model")));
              if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(base_model.get()))
                sim->initialize_simulator (this->get_simulator());

              prm.enter_subsection("Function");
              try
                {
                  function.parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Material model.Melting rate function.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
              prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        /* After parsing the parameters for depth dependent, it is essential to parse
        parameters related to the base model. */
        base_model->parse_parameters(prm);
        this->model_dependence = base_model->get_model_dependence();

        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("Material model melting rate function only "
                               "works if there is a compositional field called porosity."));
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltingRateFunction,
                                   "melting rate function",
                                   "A material model that can be derived from any of the other "
                                   "material model and that will replace the maelting rate by a "
                                   "function that can be chosen as an input parameter "
                                   "(units: 1/s or 1/yr).")

  }
}
