#include <aspect/material_model/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    Interface<dim>::~Interface ()
    {}


    template <int dim>
    double
    Interface<dim>::viscosity_derivative (const double,
                                          const double,
                                          const Point<dim> &,
                                          const NonlinearDependence::Dependence dependence) const
    {
      Assert (viscosity_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      return 0;
    }


    template <int dim>
    double
    Interface<dim>::density_derivative (const double,
                                        const double,
                                        const Point<dim> &,
                                        const NonlinearDependence::Dependence dependence) const
    {
      Assert (density_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::compressibility_derivative (const double,
                                                const double,
                                                const Point<dim> &,
                                                const NonlinearDependence::Dependence dependence) const
    {
      Assert (compressibility_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::specific_heat_derivative (const double,
                                              const double,
                                              const Point<dim> &,
                                              const NonlinearDependence::Dependence dependence) const
    {
      Assert (specific_heat_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::thermal_conductivity_derivative (const double,
                                                     const double,
                                                     const Point<dim> &,
                                                     const NonlinearDependence::Dependence dependence) const
    {
      Assert (thermal_conductivity_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      return 0;
    }


    template <int dim>
    void
    Interface<dim>::
    declare_parameters (dealii::ParameterHandler &prm)
    {}


    template <int dim>
    void
    Interface<dim>::parse_parameters (dealii::ParameterHandler &prm)
    {}


// -------------------------------- Deal with registering material models and automating
// -------------------------------- their setup and selection at run time

    namespace
    {
      internal::Plugins::PluginList<Interface<deal_II_dimension> > registered_plugins;
    }



    template <int dim>
    void
    register_material_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> *(*factory_function) ())
    {
      registered_plugins.register_plugin (name,
                                          description,
                                          declare_parameters_function,
                                          factory_function);
    }


    template <int dim>
    Interface<dim> *
    create_material_model (ParameterHandler &prm)
    {
      std::string model_name;
      prm.enter_subsection ("Material model");
      {
        model_name = prm.get ("Model name");
      }
      prm.leave_subsection ();

      return registered_plugins.create_plugin (model_name, prm);
    }


    template <int dim>
    double
    Interface<dim>::
    seismic_Vp (double dummy1,
                double dummy2) const
    {
      return -1.0;
    }


    template <int dim>
    double
    Interface<dim>::
    seismic_Vs (double dummy1,
                double dummy2) const
    {
      return -1.0;
    }


    template <int dim>
    unsigned int
    Interface<dim>::
    thermodynamic_phase (double dummy1,
                         double dummy2) const
    {
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double pressure,
                                   const Point<dim> &position) const
    {
      return (-1./density(temperature, pressure, position)
              *
              density_derivative(temperature, pressure, position, NonlinearDependence::temperature));
    }

    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the actual entry in the parameter file
      prm.enter_subsection ("Material model");
      {
        const std::string pattern_of_names
          = registered_plugins.get_pattern_of_names ();
        prm.declare_entry ("Model name", "",
                           Patterns::Selection (pattern_of_names),
                           "Select one of the following models:\n\n"
                           +
                           registered_plugins.get_description_string());
      }
      prm.leave_subsection ();

      registered_plugins.declare_parameters (prm);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<deal_II_dimension> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<deal_II_dimension> >::plugins = 0;
    }
  }

  namespace MaterialModel
  {
    template class Interface<deal_II_dimension>;

    template
    void
    register_material_model<deal_II_dimension> (const std::string &,
                                                const std::string &,
                                                void ( *) (ParameterHandler &),
                                                Interface<deal_II_dimension> *( *) ());

    template
    Interface<deal_II_dimension> *
    create_material_model<deal_II_dimension> (ParameterHandler &prm);
  }
}
