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


#ifndef _aspect_prescribed_dilation_interface_h
#define _aspect_prescribed_dilation_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/material_model/interface.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

#include <boost/core/demangle.hpp>
#include <typeinfo>


namespace aspect
{
  template <int dim> class SimulatorAccess;
  /**
   * A namespace in which we define everything that has to do with defining
   * the Prescribed dilation.
   *
   * @ingroup PrescribedDilation
   */
  namespace PrescribedDilation
  {

    template <int dim>
    class PrescribedDilationOutputs
    {
      public:
        /**
         * Constructor. Initialize the various arrays of this structure with the
         * given number of quadrature points.
         *
         * @param n_points The number of quadrature points for which output
         * quantities will be provided.
         */
        PrescribedDilationOutputs (const unsigned int n_points);

        /**
         * Copy constructor. This constructor copies all data members of the
         * source object.
         */
        PrescribedDilationOutputs (const PrescribedDilationOutputs &source) = default;

        /**
         * Move constructor. This constructor simply moves all members.
         */
        PrescribedDilationOutputs (PrescribedDilationOutputs &&)  noexcept = default;

        /**
         * Copy operator. Copying these objects is expensive, and consequently
         * prohibited.
         */
        PrescribedDilationOutputs &operator= (const PrescribedDilationOutputs &source) = delete;

        /**
         * Move operator.
         */
        PrescribedDilationOutputs &operator= (PrescribedDilationOutputs &&)  noexcept = default;

        /**
         * Function that returns the number of points at which
         * the dilation is to be evaluated.
         */
        unsigned int n_evaluation_points() const;

        /**
         * Function that set zeroes to all values in this object
         */
        void reset();

        /**
         * Vector containing the dilation at each evaluation point.
         * The physical unit of this quantity is 1/s (it is the density of volume change rate)
         */
        std::vector<double> dilation;
    };


    /**
     * A base class for parameterizations that prescribe the dilation.
     * If any Prescribed dilation plugin is selected by the user, then
     * all of the following happens:
     *
     * 1. In the continuity equation, the source term is added to the right hand side.
     *    The value we put on the RHS is computed by PrescribedDilation::evaluate function
     *    and passed to assembler in PrescribedDilationOutputs::dilation structure.
     * 2. In the momentum equation, the term grad(2/3 eta div v) is added to the left
     *    hand side. This term ensures we compute the deviatoric part of the strain rate
     *    tensor. The weak form of this term is -2/3 eta div_phi_u[i] div_phi_u[j].
     * 3. Sets do_pressure_rhs_compatibility_modification to true (if no open boundary is detected).
     *
     * If more than one prescribed dilation plugin is selected, the results are summed.
     *
     * @ingroup PrescribedDilation
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * Compute the dilation and store in PrescribedDilationOutputs structure
         */
        virtual void
        evaluate ( const aspect::MaterialModel::MaterialModelInputs<dim> &in,
                   PrescribedDilationOutputs<dim> &out) const = 0;
    };


    /**
     * A class that manages all objects that provide functionality to the
     * prescribed dilation.
     *
     * @ingroup Prescribed dilation
     */
    template <int dim>
    class Manager : public Plugins::ManagerBase<Interface<dim>>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Declare the parameters of all known dilation plugins, as
         * well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);


        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which prescribed dilation objects will be created; then
         * let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm) override;


        /**
         * Compute the dilation and store in PrescribedDilationOutputs structure
         */
        void
        evaluate ( const aspect::MaterialModel::MaterialModelInputs<dim> &in,
                   PrescribedDilationOutputs<dim> &out) const;


        /**
         * A function that is used to register dilation model objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new dilation plugin
         * class.
         *
         * @param name A string that identifies the dilation model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this dilation model wants to read
         * from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this dilation model.
         *
         * @ingroup PrescribedDilation
         */
        static
        void
        register_dilation_model (const std::string &name,
                                 const std::string &description,
                                 void (*declare_parameters_function) (ParameterHandler &),
                                 std::unique_ptr<Interface<dim>> (*factory_function) ());


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
        DeclException1 (ExcDilationModelNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered dilation model objects.");
    };


    /**
     * Return a string that consists of the names of dilation models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a dilation model, register it with the
     * aspect::PrescribedDilation::Manager class.
     *
     * @ingroup PrescribedDilation
     */
#define ASPECT_REGISTER_DILATION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_DILATION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedDilation::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::PrescribedDilation::Manager<2>::register_dilation_model, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::PrescribedDilation::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::PrescribedDilation::Manager<3>::register_dilation_model, \
                                name, description); \
  }
  }
}


#endif
