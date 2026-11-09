// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/prescribed_dilation/interface.h>

#include <deal.II/base/signaling_nan.h>

#include <tuple>
#include <list>


namespace aspect
{
  namespace PrescribedDilation
  {
    // --------------------- PrescribedDilationOutputs --------------------



    template <int dim>
    PrescribedDilationOutputs<dim>::PrescribedDilationOutputs (const unsigned int n_points) :
      dilation (n_points)
    {}

    template <int dim>
    unsigned int
    PrescribedDilationOutputs<dim>::n_evaluation_points () const
    {
      return dilation.size();
    }

    template <int dim>
    void
    PrescribedDilationOutputs<dim>::reset ()
    {
      for ( double &d: dilation )
        d = 0.0;
    }





    // ------------------------------ Manager -----------------------------



    namespace
    {
      std::tuple
      <aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::UnusablePluginList,
      aspect::internal::Plugins::PluginList<Interface<2>>,
      aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
    }


    template <int dim>
    void
    Manager<dim>::register_dilation_model (const std::string &name,
                                           const std::string &description,
                                           void (*declare_parameters_function) (ParameterHandler &),
                                           std::unique_ptr<Interface<dim>> (*factory_function) ())
    {
      std::get<dim>(registered_plugins).register_plugin (name,
                                                         description,
                                                         declare_parameters_function,
                                                         factory_function);
    }


    template <int dim>
    void
    Manager<dim>::parse_parameters (ParameterHandler &prm)
    {
      // find out which plugins are requested and the various other
      // parameters we declare here
      prm.enter_subsection ("Prescribed dilation");
      {
        this->plugin_names
          = Utilities::split_string_list(prm.get("List of model names"));

        AssertThrow(Utilities::has_unique_entries(this->plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Prescribed dilation/List of model names' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));
      }
      prm.leave_subsection ();

      // go through the list, create objects and let them parse
      // their own parameters
      for (auto &model_name : this->plugin_names)
        {
          this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                             .create_plugin (model_name,
                                                             "Prescribed dilation::Model names"));

          if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*this->plugin_objects.back()))
            sim->initialize_simulator (this->get_simulator());

          this->plugin_objects.back()->parse_parameters (prm);
          this->plugin_objects.back()->initialize ();
        }
    }


    template <int dim>
    void
    Manager<dim>::evaluate ( const aspect::MaterialModel::MaterialModelInputs<dim> &in,
                             PrescribedDilationOutputs<dim> &out) const
    {
      PrescribedDilationOutputs<dim> dilation_outputs (in.n_evaluation_points());
      out.reset();

      for (const auto &dilation_model : this->plugin_objects)
        {
          dilation_model->evaluate (in, dilation_outputs);

          for (unsigned int q = 0; q < dilation_outputs.n_evaluation_points(); ++q)
            out.dilation[q] += dilation_outputs.dilation[q];
        }
    }



    template <int dim>
    void
    Manager<dim>::declare_parameters (ParameterHandler &prm)
    {
      // declare the actual entry in the parameter file
      prm.enter_subsection ("Prescribed dilation");
      {
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of model names",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of dilation models that "
                          "will be used to calculate the source term in the continuity "
                          "equation.\n\n"
                          "The following dilation models are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());
      }
      prm.leave_subsection ();

      std::get<dim>(registered_plugins).declare_parameters (prm);
    }



    template <int dim>
    void
    Manager<dim>::write_plugin_graph (std::ostream &out)
    {
      std::get<dim>(registered_plugins).write_plugin_graph ("Prescribed dilation interface",
                                                            out);
    }



    template <int dim>
    std::string
    get_valid_model_names_pattern ()
    {
      return std::get<dim>(registered_plugins).get_pattern_of_names ();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace PrescribedDilation
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>; \
  template class PrescribedDilationOutputs<dim>; \
  \
  template \
  std::string \
  get_valid_model_names_pattern<dim> ();

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
