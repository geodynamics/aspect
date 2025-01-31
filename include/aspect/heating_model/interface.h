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


#ifndef _aspect_heating_model_interface_h
#define _aspect_heating_model_interface_h

#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/core/demangle.hpp>
#include <typeinfo>


namespace aspect
{
  template <int dim> class SimulatorAccess;
  /**
   * A namespace in which we define everything that has to do with defining
   * the heating model.
   *
   * @ingroup HeatingModels
   */
  namespace HeatingModel
  {
    /**
     * A data structure with the output field of the
     * HeatingModel::Interface::evaluate() function. The vectors are the
     * values at the different positions given by
     * MaterialModelInputs::position.
     */
    struct HeatingModelOutputs
    {
      /**
       * Constructor. Initialize the various arrays of this structure with the
       * given number of quadrature points and (finite element) components.
       *
       * @param n_points The number of quadrature points for which input
       * quantities will be provided.
       * @param n_comp The number of vector quantities (in the order in which
       * the Introspection class reports them) for which input will be
       * provided.
       */
      HeatingModelOutputs (const unsigned int n_points,
                           const unsigned int n_comp);

      /**
       * All source terms of the temperature equation that represent a
       * continuous process leading to a rate of change in temperature
       * at the given position. This includes shear heating, adiabatic
       * heating, radiogenic heat production, or any other heating rates
       * on the right hand side of the energy equation.
       */
      std::vector<double> heating_source_terms;

      /**
       * The source terms of the temperature equation that represent
       * fast changes in temperature (compared to the advection time
       * scale), for example due to reactions, at the given position.
       * On the advection time scale, these reactions might look like
       * instantaneous changes in temperatures. This includes for example
       * latent heat of melt.
       *
       * These reaction rates are only used in the operator_splitting nonlinear
       * solver scheme, which allows it to solve reactions of compositional
       * fields and temperature decoupled from the advection, and using a
       * different time step size.
       * In this case, they are used in addition to (and independent from) any
       * heating_source_terms that a heating model defines, which are assembled
       * as usual. For any other solver scheme, these values are ignored.
       *
       * In contrast to the heating source terms, these terms are actual changes
       * in temperature (units K/s or K/yr) rather than changes in energy.
       */
      std::vector<double> rates_of_temperature_change;

      /**
       * Left hand side contribution of latent heat; this is added to the
       * $\rho C_p$ term on the left hand side of the energy equation.
       */
      std::vector<double> lhs_latent_heat_terms;

      /**
       * Reset function. Resets all of the values in the heating model
       * outputs to their uninitialized values (NaN for the source and latent
       * heat terms, 0 for the rates of temperature change).
       */
      void
      reset ();
    };

    /**
     * A base class for parameterizations of heating models.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Function to compute the heating terms in @p heating_model_outputs
         * given the inputs in @p material_model_inputs and the outputs of the
         * material model in @p material_model_outputs.
         * All parts of the @p heating_model_outputs structure have to be
         * filled, heating_source_terms with the value of the heating rate and
         * lhs_latent_heat_terms with the part of the latent heat that depends
         * on the temperature change (and thus ends up on the left hand side of
         * the temperature equation) at each quadrature point as defined in
         * @p material_model_inputs, setting them to zero if they are not to
         * be used in the computation.
         */
        virtual
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const = 0;

        /**
         * Allow the heating model to attach additional material model outputs.
         * The default implementation of this function does not add any
         * outputs. Consequently, derived classes do not have to overload
         * this function if they do not need any additional outputs.
         */
        virtual
        void
        create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const;

        /**
         * Allow the heating model to attach additional material model inputs
         * it needs. The default implementation of this function does not add any
         * inputs. Consequently, derived classes do not have to overload
         * this function if they do not need any additional inputs.
         */
        virtual
        void
        create_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &inputs) const;
    };


    /**
     * A class that manages all objects that provide functionality to the
     * heating models.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Returns true if the adiabatic heating plugin is found in the
         * list of active heating models.
         */
        bool
        adiabatic_heating_enabled() const;

        /**
         * Returns true if the shear heating plugin is found in the
         * list of active heating models.
         */
        bool
        shear_heating_enabled() const;

        /**
         * Declare the parameters of all known heating plugins, as
         * well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);


        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which heating model objects will be created; then
         * let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * A function that calls the evaluate function of all the individual
         * heating models and adds up the values of the individual heating
         * model outputs.
         */
        void
        evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                  const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
                  HeatingModel::HeatingModelOutputs &heating_model_outputs) const;

        /**
         * Allow the heating model plugins to attach additional material model
         * inputs and outputs by looping over all registered plugins and calling
         * their respective member functions.
         */
        virtual
        void
        create_additional_material_model_inputs_and_outputs(MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
                                                            MaterialModel::MaterialModelOutputs<dim> &material_model_outputs) const;


        /**
         * A function that is used to register heating model objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new heating plugin
         * class.
         *
         * @param name A string that identifies the heating model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this heating model wants to read
         * from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this heating model.
         *
         * @ingroup HeatingModels
         */
        static
        void
        register_heating_model (const std::string &name,
                                const std::string &description,
                                void (*declare_parameters_function) (ParameterHandler &),
                                std::unique_ptr<Interface<dim>> (*factory_function) ());


        /**
         * Return a list of names of all heating models currently used in the
         * computation, as specified in the input file.
         *
         * @deprecated Use Plugins::ManagerBase::get_active_plugin_names()
         *   instead.
         */
        DEAL_II_DEPRECATED
        const std::vector<std::string> &
        get_active_heating_model_names () const;

        /**
         * Return a list of pointers to all heating models currently used in the
         * computation, as specified in the input file.
         *
         * @deprecated Use Plugins::ManagerBase::get_active_plugin_names()
         *   instead.
         */
        DEAL_II_DEPRECATED
        const std::list<std::unique_ptr<Interface<dim>>> &
        get_active_heating_models () const;

        /**
         * Go through the list of all heating models that have been selected
         * in the input file (and are consequently currently active) and return
         * true if one of them has the desired type specified by the template
         * argument.
         *
         * This function can only be called if the given template type (the first template
         * argument) is a class derived from the Interface class in this namespace.
         *
         * @deprecated Instead of this function, use the
         *   Plugins::ManagerBase::has_matching_active_plugin() and
         *   Plugins::ManagerBase::get_matching_active_plugin() functions of the base
         *   class of the current class.
         */
        template <typename HeatingModelType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,HeatingModelType>::value>>
        DEAL_II_DEPRECATED
        bool
        has_matching_heating_model () const;

        /**
         * Go through the list of all heating models that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be cast to that type. If so, return a reference
         * to it. If no heating model is active that matches the given type,
         * throw an exception.
         *
         * This function can only be called if the given template type (the first template
         * argument) is a class derived from the Interface class in this namespace.
         *
         * @deprecated Instead of this function, use the
         *   Plugins::ManagerBase::has_matching_active_plugin() and
         *   Plugins::ManagerBase::get_matching_active_plugin() functions of the base
         *   class of the current class.
         */
        template <typename HeatingModelType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,HeatingModelType>::value>>
        DEAL_II_DEPRECATED
        const HeatingModelType &
        get_matching_heating_model () const;


        /**
         * For the current plugin subsystem, write a connection graph of all of the
         * plugins we know about, in the format that the
         * programs dot and neato understand. This allows for a visualization of
         * how all of the plugins that ASPECT knows about are interconnected, and
         * connect to other parts of the ASPECT code.
         *
         * @param output_stream The stream to write the output to.
         */
        static
        void
        write_plugin_graph (std::ostream &output_stream);

        /**
         * Exception.
         */
        DeclException1 (ExcHeatingModelNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered heating model objects.");
    };



    template <int dim>
    template <typename HeatingModelType, typename>
    inline
    bool
    Manager<dim>::has_matching_heating_model () const
    {
      return this->template has_matching_active_plugin<HeatingModelType>();
    }


    template <int dim>
    template <typename HeatingModelType, typename>
    inline
    const HeatingModelType &
    Manager<dim>::get_matching_heating_model () const
    {
      return this->template get_matching_active_plugin<HeatingModelType>();
    }


    /**
     * Return a string that consists of the names of heating models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a heating model, register it with the
     * aspect::HeatingModel::Manager class.
     *
     * @ingroup HeatingModels
     */
#define ASPECT_REGISTER_HEATING_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_HEATING_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::HeatingModel::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::HeatingModel::Manager<2>::register_heating_model, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::HeatingModel::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::HeatingModel::Manager<3>::register_heating_model, \
                                name, description); \
  }
  }
}


#endif
