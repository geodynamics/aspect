/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/std_cxx1x/tuple.h>

#include <list>


namespace aspect
{
  namespace MaterialModel
  {
    namespace NonlinearDependence
    {
      bool
      identifies_single_variable(const Dependence dependence)
      {
        return ((dependence == temperature)
                ||
                (dependence == pressure)
                ||
                (dependence == strain_rate)
                ||
                (dependence == compositional_fields));
      }

    }

    template <int dim>
    Interface<dim>::~Interface ()
    {}

    template <int dim>
    void
    Interface<dim>::initialize ()
    {}

    template <int dim>
    void
    Interface<dim>::update ()
    {}



    template <int dim>
    double
    Interface<dim>::reference_thermal_expansion_coefficient () const
    {
      Assert(false, ExcMessage("Implement individual functions or evaluate() in material model."));
      return 1.0;
    }


    template <int dim>
    double
    Interface<dim>::viscosity_derivative (const double,
                                          const double,
                                          const std::vector<double> &, /*composition*/
                                          const Point<dim> &,
                                          const NonlinearDependence::Dependence dependence) const
    {
      Assert (viscosity_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      Assert (NonlinearDependence::identifies_single_variable(dependence) == true,
              ExcMessage ("The given dependence must identify a single variable!"));
      return 0;
    }


    template <int dim>
    double
    Interface<dim>::density_derivative (const double,
                                        const double,
                                        const std::vector<double> &, /*composition*/
                                        const Point<dim> &,
                                        const NonlinearDependence::Dependence dependence) const
    {
      Assert (density_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      Assert (NonlinearDependence::identifies_single_variable(dependence) == true,
              ExcMessage ("The given dependence must identify a single variable!"));
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::compressibility_derivative (const double,
                                                const double,
                                                const std::vector<double> &, /*composition*/
                                                const Point<dim> &,
                                                const NonlinearDependence::Dependence dependence) const
    {
      Assert (compressibility_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      Assert (NonlinearDependence::identifies_single_variable(dependence) == true,
              ExcMessage ("The given dependence must identify a single variable!"));
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::specific_heat_derivative (const double,
                                              const double,
                                              const std::vector<double> &, /*composition*/
                                              const Point<dim> &,
                                              const NonlinearDependence::Dependence dependence) const
    {
      Assert (specific_heat_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      Assert (NonlinearDependence::identifies_single_variable(dependence) == true,
              ExcMessage ("The given dependence must identify a single variable!"));
      return 0;
    }

    template <int dim>
    double
    Interface<dim>::thermal_conductivity_derivative (const double,
                                                     const double,
                                                     const std::vector<double> &, /*composition*/
                                                     const Point<dim> &,
                                                     const NonlinearDependence::Dependence dependence) const
    {
      Assert (thermal_conductivity_depends_on(dependence) == false,
              ExcMessage ("For a model declaring a certain dependence, "
                          "the partial derivatives have to be implemented."));
      Assert (NonlinearDependence::identifies_single_variable(dependence) == true,
              ExcMessage ("The given dependence must identify a single variable!"));
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
      std_cxx1x::tuple
      <void *,
      void *,
      internal::Plugins::PluginList<Interface<2> >,
      internal::Plugins::PluginList<Interface<3> > > registered_plugins;
    }



    template <int dim>
    void
    register_material_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> *(*factory_function) ())
    {
      std_cxx1x::get<dim>(registered_plugins).register_plugin (name,
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

      Interface<dim> *plugin = std_cxx1x::get<dim>(registered_plugins).create_plugin (model_name,
                                                                                      "Material model::Model name");
      return plugin;
    }


    template <int dim>
    double
    Interface<dim>::
    viscosity_ratio (const double temperature,
                     const double pressure,
                     const std::vector<double>    &compositional_fields,
                     const SymmetricTensor<2,dim> &strain_rate,
                     const Point<dim> &position) const
    {
      return 1.0;
    }


    template <int dim>
    double
    Interface<dim>::
    seismic_Vp (double dummy1,
                double dummy2,
                const std::vector<double> &, /*composition*/
                const Point<dim> &dummy3) const
    {
      return -1.0;
    }


    template <int dim>
    double
    Interface<dim>::
    seismic_Vs (double dummy1,
                double dummy2,
                const std::vector<double> &, /*composition*/
                const Point<dim> &dummy3) const
    {
      return -1.0;
    }


    template <int dim>
    unsigned int
    Interface<dim>::
    thermodynamic_phase (double dummy1,
                         double dummy2,
                         const std::vector<double> & /*composition*/) const
    {
      return 0;
    }


    template <int dim>
    void
    declare_parameters (ParameterHandler &prm)
    {
      // declare the actual entry in the parameter file
      prm.enter_subsection ("Material model");
      {
        const std::string pattern_of_names
          = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();
        try
          {
            prm.declare_entry ("Model name", "",
                               Patterns::Selection (pattern_of_names),
                               "Select one of the following models:\n\n"
                               +
                               std_cxx1x::get<dim>(registered_plugins).get_description_string());
          }
        catch (const ParameterHandler::ExcValueDoesNotMatchPattern &)
          {
            // ignore the fact that the default value for this parameter
            // does not match the pattern
          }
      }
      prm.leave_subsection ();

      std_cxx1x::get<dim>(registered_plugins).declare_parameters (prm);
    }

    template <int dim>
    Interface<dim>::MaterialModelInputs::MaterialModelInputs(unsigned int n_points, unsigned int n_comp)
    {
      position.resize(n_points);
      temperature.resize(n_points);
      pressure.resize(n_points);
      composition.resize(n_points);
      for (unsigned int i=0; i<n_points; ++i)
        composition[i].resize(n_comp);
      strain_rate.resize(n_points);
    }

    template <int dim>
    Interface<dim>::MaterialModelOutputs::MaterialModelOutputs(const unsigned int n_points,
                                                               const unsigned int n_comp)
    {
      viscosities.resize(n_points);
      densities.resize(n_points);
      thermal_expansion_coefficients.resize(n_points);
      specific_heat.resize(n_points);
      thermal_conductivities.resize(n_points);
      compressibilities.resize(n_points);
      entropy_derivative_pressure.resize(n_points);
      entropy_derivative_temperature.resize(n_points);
      reaction_terms.resize(n_points);
      for (unsigned int i=0; i<n_points; ++i)
        reaction_terms[i].resize(n_comp);
    }



    template <int dim>
    double
    InterfaceCompatibility<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double pressure,
                                   const std::vector<double> &compositional_fields,
                                   const Point<dim> &position) const
    {
      return (-1./density(temperature, pressure, compositional_fields, position)
              *
              this->density_derivative(temperature, pressure, compositional_fields, position, NonlinearDependence::temperature));
    }


    template <int dim>
    double
    InterfaceCompatibility<dim>::
    entropy_derivative (const double temperature,
                        const double pressure,
                        const std::vector<double> &compositional_fields,
                        const Point<dim> &position,
                        const NonlinearDependence::Dependence dependence) const
    {
      return 0.0;
    }


    template <int dim>
    double
    InterfaceCompatibility<dim>::
    reaction_term (const double temperature,
                   const double pressure,
                   const std::vector<double> &compositional_fields,
                   const Point<dim> &position,
                   const unsigned int compositional_variable) const
    {
      return 0.0;
    }


    template <int dim>
    void
    InterfaceCompatibility<dim>::evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                                          typename Interface<dim>::MaterialModelOutputs &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
        {
          // as documented, if the strain rate array is empty, then do not compute the
          // viscosities
          if (in.strain_rate.size() > 0)
            out.viscosities[i]                  = viscosity                     (in.temperature[i], in.pressure[i], in.composition[i], in.strain_rate[i], in.position[i]);

          out.densities[i]                      = density                       (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.thermal_expansion_coefficients[i] = thermal_expansion_coefficient (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.specific_heat[i]                  = specific_heat                 (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.thermal_conductivities[i]         = thermal_conductivity          (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.compressibilities[i]              = compressibility               (in.temperature[i], in.pressure[i], in.composition[i], in.position[i]);
          out.entropy_derivative_pressure[i]    = entropy_derivative            (in.temperature[i], in.pressure[i], in.composition[i], in.position[i], NonlinearDependence::pressure);
          out.entropy_derivative_temperature[i] = entropy_derivative            (in.temperature[i], in.pressure[i], in.composition[i], in.position[i], NonlinearDependence::temperature);
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            out.reaction_terms[i][c]            = reaction_term                 (in.temperature[i], in.pressure[i], in.composition[i], in.position[i], c);
        }
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
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<2> >::plugins = 0;

      template <>
      std::list<internal::Plugins::PluginList<MaterialModel::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<MaterialModel::Interface<3> >::plugins = 0;
    }
  }

  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template class InterfaceCompatibility<dim>; \
  \
  template \
  void \
  register_material_model<dim> (const std::string &, \
                                const std::string &, \
                                void ( *) (ParameterHandler &), \
                                Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  Interface<dim> * \
  create_material_model<dim> (ParameterHandler &prm);

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
