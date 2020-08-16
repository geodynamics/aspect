/*
  Copyright (C) 2011 - 2020 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_interface_h
#define _aspect_material_model_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>
#include <aspect/material_model/utilities.h>

#include <deal.II/base/point.h>

// Work around an incorrect instantiation in qprojector.h of deal.II 9.2.0,
// which requires including qprojector.h before quadrature.h (and not
// after). This file doesn't actually need qprojector.h, so the include can be
// removed when we require 9.3.. For more info see
// https://github.com/geodynamics/aspect/issues/3728
#if !DEAL_II_VERSION_GTE(9,3,0)
#include <deal.II/base/qprojector.h>
#endif
#include <deal.II/base/quadrature.h>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/base/signaling_nan.h>

namespace aspect
{
  template <int dim>
  struct Introspection;


  /**
   * A namespace in which we define everything that has to do with modeling
   * convecting material, including descriptions of material parameters such
   * as viscosities, densities, etc.
   *
   * @ingroup MaterialModels
   */
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * A namespace whose enum members are used in querying the nonlinear
     * dependence of physical parameters on other solution variables.
     */
    namespace NonlinearDependence
    {
      /**
       *
       * An enum whose members are used in querying the nonlinear dependence
       * of physical parameters on other solution variables.
       *
       * The values of this enum are used in the
       * NonlinearDependence::model_dependence, queried by get_model_dependence()
       * to see if a coefficient like the viscosity, depends on the
       * temperature, pressure, strain rate, or compositional field value.
       * Because the values of the enum are chosen so that they represent
       * single bits in an integer, the result here is a number that can be
       * represented in base-2 as 101 (the number 100=4 for the strain rate and
       * 001=1 for the temperature).
       *
       * To query nonlinear dependence of a coefficient on any other variable,
       * you can use
       * @code
       *   material_model.get_model_dependence();
       * @endcode
       * and compare the result to NonlinearDependence::Dependence::Variable.
       */
      enum Dependence
      {
        uninitialized        = 0,

        none                 = 1,
        temperature          = 2,
        pressure             = 4,
        strain_rate          = 8,
        compositional_fields = 16,

        any_variable         = temperature | pressure | strain_rate | compositional_fields
      };


      /**
       * Provide an operator that or's two Dependence variables.
       */
      inline Dependence operator | (const Dependence d1,
                                    const Dependence d2)
      {
        return Dependence(static_cast<int>(d1) | static_cast<int>(d2));
      }

      inline Dependence operator |= (Dependence &d1,
                                     const Dependence d2)
      {
        d1 = (d1 | d2);
        return d1;
      }

      /**
       * A structure that, for every output variable of a material model,
       * describes which input variable it depends on.
       */
      struct ModelDependence
      {
        /**
         * A field that describes which input variable the viscosity
         * of a material model depends on.
         */
        Dependence viscosity;

        /**
         * A field that describes which input variable the density
         * of a material model depends on.
         */
        Dependence density;

        /**
         * A field that describes which input variable the compressibility
         * of a material model depends on.
         */
        Dependence compressibility;

        /**
         * A field that describes which input variable the specific heat
         * of a material model depends on.
         */
        Dependence specific_heat;

        /**
         * A field that describes which input variable the thermal conductivity
         * of a material model depends on.
         */
        Dependence thermal_conductivity;

        /**
         * Default constructor. Sets all dependencies to invalid values in
         * order to ensure that material models really correctly specify
         * which input variables their output variables depend on.
         */
        ModelDependence ();
      };

      /**
       * Return whether the given argument @p dependence identifies a single
       * variable (e.g., the pressure, the temperature, etc) or a combination
       * of variables. Technically, this corresponds to the question of
       * whether there is exactly one bit set in the argument.
       *
       * @return true if yes, false otherwise.
       */
      bool
      identifies_single_variable(const Dependence dependence);
    }

    /**
    * A namespace whose enum members are used in querying which material
    * properties should be computed.
    */
    namespace MaterialProperties
    {
      /**
       * An enum whose members identify material model output
       * properties.
       *
       * Because the values of the enum are chosen so that they represent
       * single bits in an integer, the result here is a number that can be
       * represented in base-2 as 110 (the number 100=4 for the density and
       * 010=2 for the viscosity).
       */
      enum Property
      {
        uninitialized                  = 0,

        none                           = 1,
        viscosity                      = 2,
        density                        = 4,
        thermal_expansion_coefficient  = 8,
        specific_heat                  = 16,
        thermal_conductivity           = 32,
        compressibility                = 64,
        entropy_derivative_pressure    = 128,
        entropy_derivative_temperature = 256,
        reaction_terms                 = 512,

        equation_of_state_properties   = density |
                                         thermal_expansion_coefficient |
                                         specific_heat |
                                         compressibility |
                                         entropy_derivative_pressure |
                                         entropy_derivative_temperature,
        all_properties                 = equation_of_state_properties |
                                         viscosity |
                                         reaction_terms
      };

      /**
       * Provide an operator that or's two Property variables. This allows to
       * combine more than one property in a single variable.
       */
      inline Property operator | (const Property d1,
                                  const Property d2)
      {
        return Property(static_cast<int>(d1) | static_cast<int>(d2));
      }

      inline Property operator |= (Property &d1,
                                   const Property d2)
      {
        return (d1 | d2);
      }
    }

    template <int dim> class AdditionalMaterialInputs;

    /**
     * A data structure with all inputs for the
     * MaterialModel::Interface::evaluate() method. The vectors all have the
     * same length and refer to different evaluation points (given in
     * #position).
     */
    template <int dim>
    struct MaterialModelInputs
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
      MaterialModelInputs(const unsigned int n_points,
                          const unsigned int n_comp);

      /**
       * Constructor. Initialize the arrays of the structure with the number
       * of points in the `input_data` structure, and fills them appropriately.
       *
       * @param input_data The data used to populate the material model input quantities.
       * @param introspection A reference to the simulator introspection object.
       * @param use_strain_rate Whether to compute the strain rates.
       */
      MaterialModelInputs(const DataPostprocessorInputs::Vector<dim> &input_data,
                          const Introspection<dim> &introspection,
                          const bool use_strain_rate = true);


      /**
       * Constructor. Initializes the various arrays of this
       * structure with the FEValues and introspection objects and
       * the solution_vector. This constructor calls the function
       * reinit to populate the newly created arrays.
       *
       * @param fe_values An FEValuesBase object used to evaluate the finite elements.
       * @param cell The currently active cell for the fe_values object.
       * @param introspection A reference to the simulator introspection object.
       * @param solution_vector The finite element vector from which to construct the inputs.
       * @param use_strain_rates Whether to compute the strain rates.
       */
      MaterialModelInputs(const FEValuesBase<dim,dim> &fe_values,
                          const typename DoFHandler<dim>::active_cell_iterator &cell,
                          const Introspection<dim> &introspection,
                          const LinearAlgebra::BlockVector &solution_vector,
                          const bool use_strain_rates = true);

      /**
       * Copy constructor. This constructor copies all data members of the
       * source object except for the additional input data (of type
       * AdditionalMaterialInputs) pointers, stored in the
       * `source.additional_inputs` member variable.
       *
       * This is because these pointers can not be copied (they
       * are unique to the @p source object). Since they can also not
       * be recreated without the original code that created these objects
       * in the first place, this constructor throws an exception if the
       * @p source object had any additional input data objects
       * associated with it.
       */
      MaterialModelInputs (const MaterialModelInputs &source);

      /**
       * Move constructor. This constructor simply moves all members.
       */
      MaterialModelInputs (MaterialModelInputs &&) = default;

      /**
       * Copy operator. Copying these objects is expensive and
       * consequently prohibited
       */
      MaterialModelInputs &operator= (const MaterialModelInputs &source) = delete;

      /**
       * Move operator.
       */
      MaterialModelInputs &operator= (MaterialModelInputs &&) = default;

      /**
       * Function to re-initialize and populate the pre-existing arrays
       * created by the constructor MaterialModelInputs.
       */
      void reinit(const FEValuesBase<dim,dim> &fe_values,
                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                  const Introspection<dim> &introspection,
                  const LinearAlgebra::BlockVector &solution_vector,
                  const bool use_strain_rates = true);

      /**
       * Function that returns the number of points at which
       * the material model is to be evaluated.
       */
      unsigned int n_evaluation_points() const;

      /**
       * Function that returns if the caller requests an evaluation
       * of the handed over @p property. This is optional, because calculating
       * some properties can be more expensive than the other material
       * model properties and not all are needed for all applications.
       */
      bool requests_property(const MaterialProperties::Property &property) const;

      /**
       * Vector with global positions where the material has to be evaluated
       * in evaluate().
       */
      std::vector<Point<dim> > position;

      /**
       * Temperature values at the points given in the #position vector.
       */
      std::vector<double> temperature;

      /**
       * Pressure values at the points given in the #position vector.
       */
      std::vector<double> pressure;

      /**
       * Pressure gradients at the points given in the #position vector.
       * This is important for the heating models.
       */
      std::vector<Tensor<1,dim> > pressure_gradient;

      /**
       * Velocity values at the points given in the #position vector.
       * This value is mostly important in the case of determining
       * whether material crossed a certain region (e.g. a phase boundary).
       * The timestep that is needed for this check can be requested from
       * SimulatorAccess.
       */
      std::vector<Tensor<1,dim> > velocity;

      /**
       * Values of the compositional fields at the points given in the
       * #position vector: composition[i][c] is the compositional field c at
       * point i.
       */
      std::vector<std::vector<double> > composition;

      /**
       * Strain rate at the points given in the #position vector. Only the
       * viscosity may depend on these values. This std::vector can be set to
       * size 0 if the viscosity is not needed.
       *
       * @note The strain rate is computed as $\varepsilon(\mathbf u)=\frac 12
       * (\nabla \mathbf u + \nabla \mathbf u^T)$, regardless of whether the
       * model is compressible or not. This is relevant since in some other
       * contexts, the strain rate in the compressible case is computed as
       * $\varepsilon(\mathbf u)=\frac 12 (\nabla \mathbf u + \nabla \mathbf
       * u^T) - \frac 13 \nabla \cdot \mathbf u \mathbf 1$.
       */
      std::vector<SymmetricTensor<2,dim> > strain_rate;

      /**
       * Optional cell object that contains these quadrature
       * points. This allows for evaluating properties at the cell vertices
       * and interpolating to the quadrature points, or to query the cell for
       * material ids, neighbors, or other information that is not available
       * solely from the locations. Note that not all calling functions will
       * set this cell iterator. In these cases it will be an invalid iterator
       * constructed using the default constructor, so make sure that your
       * material model either fails
       * with a proper error message, or provides an alternative calculation for
       * these cases. You can detect this with
       * @code
       * if (in.current_cell.state() == IteratorState::valid)
       * @endcode
       */
      typename DoFHandler<dim>::active_cell_iterator current_cell;

      /**
       * A member variable that stores which properties the material model
       * should compute. You can check specific properties using
       * the requests_property function and usually do not need to access
       * this variable directly. For documentation on the internal storage
       * of this variable see the documentation for MaterialProperties::Property.
       */
      MaterialProperties::Property requested_properties;

      /**
       * Vector of shared pointers to additional material model input
       * objects that can be added to MaterialModelInputs. By default,
       * no inputs are added.
       */
      std::vector<std::unique_ptr<AdditionalMaterialInputs<dim> > > additional_inputs;

      /**
       * Given an additional material model input class as explicitly specified
       * template argument, returns a pointer to this additional material model
       * input object if it is used in the current simulation.
       * If the output does not exist, a null pointer is returned.
       */
      template <class AdditionalInputType>
      AdditionalInputType *get_additional_input();

      /**
       * Constant version of get_additional_input() returning a const pointer.
       */
      template <class AdditionalInputType>
      const AdditionalInputType *get_additional_input() const;
    };


    template <int dim>     class AdditionalMaterialOutputs;


    /**
     * A data structure with the output field of the
     * MaterialModel::Interface::evaluate() function. The vectors are the
     * values at the different positions given by
     * MaterialModelInputs::position.
     */
    template <int dim>
    struct MaterialModelOutputs
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
      MaterialModelOutputs (const unsigned int n_points,
                            const unsigned int n_comp);

      /**
       * Copy constructor. This constructor copies all data members of the
       * source object except for the additional output data (of type
       * AdditionalMaterialOutputs) pointers, stored in the
       * `source.additional_outputs` member variable.
       *
       * This is because these pointers can not be copied (they
       * are unique to the @p source object). Since they can also not
       * be recreated without the original code that created these objects
       * in the first place, this constructor throws an exception if the
       * @p source object had any additional output data objects
       * associated with it.
       */
      MaterialModelOutputs (const MaterialModelOutputs &source);

      /**
       * Move constructor. This constructor simply moves all members.
       */
      MaterialModelOutputs (MaterialModelOutputs &&) = default;

      /**
       * Copy operator. Copying these objects is expensive, and consequently
       * prohibited.
       */
      MaterialModelOutputs &operator= (const MaterialModelOutputs &source) = delete;

      /**
       * Move operator.
       */
      MaterialModelOutputs &operator= (MaterialModelOutputs &&) = default;

      /**
      * Function that returns the number of points at which
      * the material model is to be evaluated.
      */
      unsigned int n_evaluation_points() const;

      /**
       * Viscosity $\eta$ values at the given positions.
       */
      std::vector<double> viscosities;

      /**
       * Density values at the given positions.
       */
      std::vector<double> densities;

      /**
       * Thermal expansion coefficients at the given positions. It is defined
       * as $\alpha = - \frac{1}{\rho} \frac{\partial\rho}{\partial T}$
       */
      std::vector<double> thermal_expansion_coefficients;

      /**
       * Specific heat at the given positions.
       */
      std::vector<double> specific_heat;

      /**
       * Thermal conductivity at the given positions.
       */
      std::vector<double> thermal_conductivities;

      /**
       * Compressibility at the given positions. The compressibility is defined
       * as $\kappa = \frac{1}{\rho} \frac{\partial\rho}{\partial p}$.
       */
      std::vector<double> compressibilities;

      /**
       * The product of the change of entropy $\Delta S$ at a phase transition
       * and the derivative of the phase function $X=X(p,T,\mathfrak c,\mathbf
       * x)$ with regard to pressure at the given positions.
       */
      std::vector<double> entropy_derivative_pressure;

      /**
       * The product of (minus) the change of entropy $-\Delta S$ at a phase
       * transition and the derivative of the phase function
       * $X=X(p,T,\mathfrak c,\mathbf x)$ with regard to temperature at the
       * given positions.
       */
      std::vector<double> entropy_derivative_temperature;

      /**
       * Change in composition due to chemical reactions at the given
       * positions. The term reaction_terms[i][c] is the change in
       * compositional field c at point i.
       *
       * The mental model behind prescribing actual changes in composition
       * rather than reaction rates is that we assume that there is always an
       * equilibrium between the compositional fields (because the time scale
       * of reactions is normally much shorter than that of convection), so
       * the quantity returned by this function is an actual change in the
       * amount of material, which is added to or subtracted from the current
       * value of the compositional field, and NOT a reaction rate. The idea
       * is, that in dependence of temperature, pressure, position and the
       * compositional fields themselves an equilibrium can be calculated, and
       * the difference between the current value and the equilibrium can be
       * added to the respective compositional field.
       *
       * For mass conservation it should ALWAYS be checked that what is
       * subtracted from one field is added to another field (and the other
       * way round) and that one never subtracts more than the actual value of
       * a field (so it does not get negative).
       *
       * This function has a default implementation that sets the reaction
       * term to zero (assuming no reactions).
       *
       * @note In cases where one has slow chemical reactions (or cases where
       * compositional fields are used to track quantities different than
       * actual compositions, for example accumulated strains in damage
       * models), models are formulated as differential equations with right
       * hand sides, not as instantaneous equations. In such cases, the
       * reaction terms (i.e., the incremental additions to the previous
       * state) are usually of the form reaction rate times time step size. To
       * implement something like this, derive your material model from
       * SimulatorAccess so you can query the time step used by the simulator
       * in order to compute the reaction increment.
       */
      std::vector<std::vector<double> > reaction_terms;

      /**
       * Vector of shared pointers to additional material model output
       * objects that can then be added to MaterialModelOutputs. By default,
       * no outputs are added.
       */
      std::vector<std::unique_ptr<AdditionalMaterialOutputs<dim> > > additional_outputs;

      /**
       * Given an additional material model output class as explicitly specified
       * template argument, returns a pointer to this additional material model
       * output object if it is used in the current simulation.
       * The output can then be filled in the MaterialModels::Interface::evaluate()
       * function. If the output does not exist, a null pointer is returned.
       */
      template <class AdditionalOutputType>
      AdditionalOutputType *get_additional_output();

      /**
       * Constant version of get_additional_output() returning a const pointer.
       */
      template <class AdditionalOutputType>
      const AdditionalOutputType *get_additional_output() const;

      /**
       * Steal the additional outputs from @p other. The destination (@p
       * this), is expected to currently have no additional outputs.
       */
      void move_additional_outputs_from(MaterialModelOutputs<dim> &other);
    };



    /**
     * A namespace in which we define how material model outputs should be
     * averaged on each cell.
     *
     * Material models compute output quantities such as the viscosity, the
     * density, etc, based on pressures, temperatures, composition, and
     * location at every quadrature point. For some models, these values vary
     * drastically from quadrature point to quadrature point, and this creates
     * difficulties both for the stability of the discretization as well as
     * for the linear solvers. Some of this can be ameliorated by averaging
     * values on every cell, although this of course reduces the ideal
     * convergence order. This namespace defines the means to achieve such
     * averaging.
     */
    namespace MaterialAveraging
    {
      /**
       * An enum to define what kind of averaging operations are implemented.
       * These are:
       *
       * - No averaging, i.e., leave the values as they were provided by the
       * material model.
       *
       * - Arithmetic averaging: Set the values of each output quantity at
       * every quadrature point to \f[ \bar x = \frac 1Q \sum_{q=1}^Q x_q \f]
       * where $x_q$ are the values at the $Q$ quadrature points.
       *
       * - Harmonic averaging: Set the values of each output quantity at every
       * quadrature point to \f[ \bar x = \left(\frac 1Q \sum_{q=1}^Q
       * \frac{1}{x_q}\right)^{-1} \f] where $x_q$ are the values at the $Q$
       * quadrature points.
       *
       * - Geometric averaging: Set the values of each output quantity at
       * every quadrature point to \f[ \bar x = \left(\prod_{q=1}^Q
       * x_q\right)^{1/Q} \f] where $x_q$ are the values at the $Q$ quadrature
       * points.
       *
       * - Pick largest: Set the values of each output quantity at every
       * quadrature point to \f[ \bar x = \max_{1\le q\le Q} x_q \f] where
       * $x_q$ are the values at the $Q$ quadrature points.
       *
       * - Project to $Q_1$: This operation takes the values at the quadrature
       * points and computes the best bi- or trilinear approximation for them.
       * In other words, it projects the values into the $Q_1$ finite element
       * space. It then re-evaluate this projection at the quadrature points.
       *
       * - Log average: Set the values of each output quantity at every
       * quadrature point to \f[ \bar x = {10}^{\frac 1Q \sum_{q=1}^Q \log_{10} x_q} \f]
       * where $x_q$ are the values at the $Q$ quadrature points.
       *
       * - Harmonic average only viscosity and project to Q1 only viscosity: Like
       * harmonic averaging and project to Q1, but only
       * applied to the viscosity.
       */
      enum AveragingOperation
      {
        none,
        arithmetic_average,
        harmonic_average,
        geometric_average,
        pick_largest,
        project_to_Q1,
        log_average,
        harmonic_average_only_viscosity,
        project_to_Q1_only_viscosity
      };


      /**
       * Return a string that represents the various averaging options laid
       * out above and that can be used in the declaration of an input
       * parameter. The options are separated by "|" so that they can be used
       * in a dealii::Patterns::Selection argument.
       */
      std::string get_averaging_operation_names ();

      /**
       * Parse a string representing one of the options returned by
       * get_averaging_operation_names(), and return the corresponding
       * AveragingOperation value.
       */
      AveragingOperation parse_averaging_operation_name (const std::string &s);

      /**
       * Given the averaging @p operation, a description of where the
       * quadrature points are located on the given cell, and a mapping,
       * perform this operation on all elements of the @p values structure.
       */
      template <int dim>
      void average (const AveragingOperation operation,
                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                    const Quadrature<dim>         &quadrature_formula,
                    const Mapping<dim>            &mapping,
                    MaterialModelOutputs<dim>     &values_out);

      /**
       * Do the requested averaging operation for one array. The
       * projection matrix argument is only used if the operation
       * chosen is project_to_Q1
       */
      void average_property (const AveragingOperation  operation,
                             const FullMatrix<double>      &projection_matrix,
                             const FullMatrix<double>      &expansion_matrix,
                             std::vector<double>           &values_out);
    }


    /**
     * Some material and heating models need more than just the basic material
     * model inputs defined in the MaterialModel::MaterialModelInputs
     * class. These additions are either for more complicated physics
     * than the basic flow model usually solved by ASPECT (for example
     * to support the melt migration functionality), or other derived
     * quantities.
     *
     * Rather than litter the MaterialModelInputs class with additional
     * fields that are not universally used, we use a mechanism by
     * which MaterialModelInputs can store a set of pointers to
     * "additional" input objects that store information such as
     * mentioned above. These pointers are all to objects whose types
     * are derived from the current class.
     *
     * The format of the additional quantities defined in derived classes
     * should be the same as for MaterialModel::MaterialModelInputs.
     */
    template <int dim>
    class AdditionalMaterialInputs
    {
      public:
        virtual ~AdditionalMaterialInputs()
        {}

        /**
         * Fill the additional inputs. Each additional input
         * has to implement their own version of this function.
         */
        virtual void
        fill (const LinearAlgebra::BlockVector &solution,
              const FEValuesBase<dim>          &fe_values,
              const Introspection<dim>         &introspection) = 0;
    };


    /**
     * Some material models can compute more than just the basic material
     * coefficients defined in the MaterialModel::MaterialModelOutputs
     * class. These additions are either for more complicated physics
     * than the basic flow model usually solved by ASPECT (for example
     * to support the melt migration functionality), or other derived
     * quantities that are not coefficients in any of the equations
     * ASPECT solves but that may be of interest for visualization
     * (for example seismic velocities).
     *
     * Rather than litter the MaterialModelOutputs class with additional
     * fields that are not universally used, we use a mechanism by
     * which MaterialModelOutputs can store a set of pointers to
     * "additional" output objects that store information such as
     * mentioned above. These pointers are all to objects whose types
     * are derived from the current class.
     *
     * If an implementation of the MaterialModel::Interface::evaluate()
     * in a class derived from MaterialModel::Interface encounters
     * a MaterialModelOutputs object that has these pointers set
     * (and if it recognizes the type of the object pointed to),
     * it should fill this set of additional output quantities.
     *
     * The format of the additional quantities defined in derived classes
     * should be the same as for MaterialModel::MaterialModelOutputs.
     */
    template <int dim>
    class AdditionalMaterialOutputs
    {
      public:
        /**
         * Destructor. Made virtual to enable storing pointers to this
         * base class.
         */
        virtual ~AdditionalMaterialOutputs() = default;

        virtual void average (const MaterialAveraging::AveragingOperation /*operation*/,
                              const FullMatrix<double>  &/*projection_matrix*/,
                              const FullMatrix<double>  &/*expansion_matrix*/)
        {}
    };


    /**
     * Some material models can compute things that are not used anywhere
     * in the physics modules of ASPECT, but that may be of interest for
     * visualization purposes. An example would be a material model that can
     * compute seismic velocities -- these are irrelevant to the rest of
     * ASPECT, but would be nice to have for postprocessing.
     *
     * This class is a base class for material models to provide this kind
     * of information. It follows the scheme laid out by
     * AdditionalMaterialModelOutputs but also provides an interface by which
     * consumers of these objects (e.g., the
     * Postprocess::Visualization::NamedAdditionalOutputs class) can query the
     * names and values material models have put into these additional
     * outputs. (Because every material model can decide by itself which --
     * if any -- additional outputs it produces, there are no standard
     * names. Consequently, the material models have to describe what
     * values and how many values they can produce.)
     *
     * This class is then this base class for additional named material model outputs
     * to be added to the MaterialModel::MaterialModelOutputs structure.
     */
    template <int dim>
    class NamedAdditionalMaterialOutputs : public AdditionalMaterialOutputs<dim>
    {
      public:
        /**
         * Base constructor.
         *
         * @param output_names A list of names for the additional output variables
         *   this object will store. The length of the list also indicates
         *   how many additional output variables objects of derived classes
         *   will store.
         */
        NamedAdditionalMaterialOutputs(const std::vector<std::string> &output_names);

        /**
         * Constructor for case where outputs are stored for a number of points.
         *
         * @param output_names A list of names for the additional output variables
         *   this object will store. The length of the list also indicates
         *   how many additional output variables objects of derived classes
         *   will store.
         * @param n_points The number of points for which to store each of the
         *   output variables.
         */
        NamedAdditionalMaterialOutputs(const std::vector<std::string> &output_names,
                                       const unsigned int n_points);

        /**
         * Destructor.
         */
        ~NamedAdditionalMaterialOutputs() override;

        /**
         * Return a reference to the vector of names of the additional
         * outputs.
         */
        const std::vector<std::string> &get_names() const;

        /**
         * Given an index as input argument, return a reference the to vector of
         * values of the additional output with that index.
         */
        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        void average (const MaterialAveraging::AveragingOperation /*operation*/,
                      const FullMatrix<double>  &/*projection_matrix*/,
                      const FullMatrix<double>  &/*expansion_matrix*/) override
        {}


        /**
         * Values for the outputs at a set of evaluation points
         * output_values[i][j] is the value of output i at point j.
         */
        std::vector<std::vector<double> > output_values;

      private:
        const std::vector<std::string> names;
    };


    /**
     * Additional output fields for the seismic velocities to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class SeismicAdditionalOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        SeismicAdditionalOutputs(const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Seismic s-wave velocities at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> vs;

        /**
         * Seismic p-wave velocities at the evaluation points passed to
         * the instance of MaterialModel::Interface::evaluate() that fills
         * the current object.
         */
        std::vector<double> vp;
    };


    /**
     * Additional output fields for reaction rates to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     *
     * These reaction rates are only used if the "operator splitting" solver scheme
     * option is enabled, which decouples the reactions between compositional
     * fields from the advection, so that different time step sizes can be used.
     * In this case, the reaction rates are used in addition to (and independent
     * from) any reaction_terms that a material model defines, which are assembled
     * as usual.
     * By default, the reaction rates are initialized with quiet_NaNs, and if
     * "operator splitting" is not enabled, these values are not used, and they
     * are expected to either remain at that value, or to not be created at all.
     *
     * In contrast to the reaction_terms, which are actual changes in composition
     * rather than reaction rates, and assume equilibrium between the compositional
     * fields, the reaction_rates defined here allow for reaction processes that
     * happen on shorter time scales than the advection, and disequilibrium reactions.
     */
    template <int dim>
    class ReactionRateOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        ReactionRateOutputs (const unsigned int n_points,
                             const unsigned int n_comp);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Reaction rates for all compositional fields at the evaluation points
         * that are passed to the instance of MaterialModel::Interface::evaluate()
         * that fills the current object.
         * reaction_rates[q][c] is the reaction rate at the evaluation point q
         * for the compositional field with the index c.
         */
        std::vector<std::vector<double> > reaction_rates;
    };


    /**
     * Additional output fields for prescribed field outputs to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     *
     * This structure is used if for one of the compositional fields employed
     * by a simulation, the advection scheme "prescribed field" is selected.
     * (See Parameters::AdvectionFieldMethod for more information.)
     * Then, while updating the compositional field, a structure of this
     * type is created, given to the material model, and the material model
     * outputs will finally be interpolated onto the corresponding
     * compositional field.
     *
     * @note This structure always has as many prescribed field
     * outputs as there are compositional fields, even if not all of them
     * are using the "prescribed field" method. It is the responsibility
     * of the individual material models to fill the correct entries.
     */
    template <int dim>
    class PrescribedFieldOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        PrescribedFieldOutputs (const unsigned int n_points,
                                const unsigned int n_comp);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Prescribed field outputs for all compositional fields at the evaluation points
         * that are passed to the instance of MaterialModel::Interface::evaluate()
         * that fills the current object.
         * prescribed_field_outputs[q][c] is the prescribed field output at the evaluation point q
         * for the compositional field with the index c.
         */
        std::vector<std::vector<double> > prescribed_field_outputs;
    };

    /**
     * Additional output fields for prescribed temperature outputs to be added to
     * the MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     *
     * This structure is used if the temperature field is solved using the
     * advection scheme "prescribed field".
     * (See Parameters::AdvectionFieldMethod for more information.)
     * Then, while updating the temperature field, a structure of this
     * type is created, given to the material model, and the material model
     * outputs will finally be interpolated onto the
     * temperature field.
     */
    template <int dim>
    class PrescribedTemperatureOutputs : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        PrescribedTemperatureOutputs (const unsigned int n_points);

        std::vector<double> get_nth_output(const unsigned int idx) const override;

        /**
         * Prescribed field outputs for the temperature field at the evaluation points
         * that are passed to the instance of MaterialModel::Interface::evaluate()
         * that fills the current object.
         * prescribed_temperature_outputs[q] is the prescribed field output at the evaluation point q
         * for the temperature.
         */
        std::vector<double> prescribed_temperature_outputs;
    };

    /**
     * A class for additional output fields to be added to the RHS of the
     * Stokes system, which can be attached to the
     * MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class AdditionalMaterialOutputsStokesRHS: public AdditionalMaterialOutputs<dim>
    {
      public:
        AdditionalMaterialOutputsStokesRHS(const unsigned int n_points)
          : rhs_u(n_points), rhs_p(n_points), rhs_melt_pc(n_points)
        {}

        ~AdditionalMaterialOutputsStokesRHS() override
        {}

        void average (const MaterialAveraging::AveragingOperation /*operation*/,
                      const FullMatrix<double>  &/*projection_matrix*/,
                      const FullMatrix<double>  &/*expansion_matrix*/) override
        {
          // TODO: not implemented
        }

        /**
         * Force tensor on the right-hand side for the conservation of
         * momentum equation (first part of the Stokes equation) in each
         * quadrature point.
         */
        std::vector<Tensor<1,dim> > rhs_u;

        /**
         * Force value for the conservation of mass equation (second Stokes
         * equation) in each quadrature point.
         */
        std::vector<double> rhs_p;

        /**
         * Force for the compaction pressure equation (when using melt
         * transport) in each quadrature point.
         */
        std::vector<double> rhs_melt_pc;
    };



    /**
     * An AdditionalOutput that allows prescribing a dilation applied to the
     * Stokes solution.
     *
     * This is typically used in a MaterialModel to add dilation when plastic
     * failure occurs as motivated by ChoiPeterson2015. If this output
     * (denoted by R below) is present and enable_prescribed_dilation==true
     * the following terms will be assembled:
     *
     * 1) $\int - (R,q)$ to the conservation of mass equation, creating
     *    $-(div u,q) = -(R,q)$.
     * 2) $\int - 2.0 / 3.0 * eta * (R, div v)$ to the RHS of the momentum
     *    equation (if the model is incompressible), otherwise this term is
     *    already present on the left side.
     */
    template <int dim>
    class PrescribedPlasticDilation : public NamedAdditionalMaterialOutputs<dim>
    {
      public:
        /**
         * Constructor
         */
        explicit PrescribedPlasticDilation (const unsigned int n_points);

        /**
         * Function for NamedAdditionalMaterialOutputs interface
         */
        virtual std::vector<double> get_nth_output(const unsigned int idx) const;

        /**
         * A scalar value per evaluation point that specifies the prescribed dilation
         * in that point.
         */
        std::vector<double> dilation;
    };



    /**
     * A class for an elastic force term to be added to the RHS of the
     * Stokes system, which can be attached to the
     * MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function.
     */
    template <int dim>
    class ElasticOutputs: public AdditionalMaterialOutputs<dim>
    {
      public:
        ElasticOutputs(const unsigned int n_points)
          : elastic_force(n_points, numbers::signaling_nan<Tensor<2,dim> >() )
        {}

        ~ElasticOutputs() override
        {}

        void average (const MaterialAveraging::AveragingOperation operation,
                      const FullMatrix<double>  &/*projection_matrix*/,
                      const FullMatrix<double>  &/*expansion_matrix*/) override
        {
          AssertThrow(operation == MaterialAveraging::AveragingOperation::none,ExcNotImplemented());
          return;
        }

        /**
         * Force tensor (elastic terms) on the right-hand side for the conservation of
         * momentum equation (first part of the Stokes equation) in each
         * quadrature point.
         */
        std::vector<Tensor<2,dim> > elastic_force;
    };



    /**
     * A base class for parameterizations of material models. Classes derived
     * from this class will need to implement functions that provide material
     * parameters such as the viscosity, density, etc, typically as a function
     * of position, temperature and pressure at that location.
     *
     * Implementing a material model requires you to override evaluate() and fill the output
     * argument struct instead of implementing the functions viscosity(),
     * density(), etc.. In this case, all other functions are being ignored.
     *
     * In all cases, model_dependence values, is_compressible(), reference_viscosity()
     * need to be implemented.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * A typedef to import the MaterialModelInputs name into the current
         * class. This typedef primarily exists as a backward compatibility
         * measure given that the referenced structure used to be a member of
         * the current class.
         */
        using MaterialModelInputs = MaterialModel::MaterialModelInputs<dim>;
        /**
         * A typedef to import the MaterialModelOutputs name into the current
         * class. This typedef primarily exists as a backward compatibility
         * measure given that the referenced structure used to be a member of
         * the current class.
         */
        using MaterialModelOutputs = MaterialModel::MaterialModelOutputs<dim>;

        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * Called at the beginning of each time step and allows the material
         * model to update internal data structures.
         */
        virtual void update ();

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */

        /**
         * Return a structure that describes how each of the model's
         * output variables (such as viscosity, density, etc) depend
         * on the input variables pressure, temperature, strain rate,
         * and compositional fields.
         */
        const NonlinearDependence::ModelDependence &
        get_model_dependence () const;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the
         * continuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        virtual bool is_compressible () const = 0;
        /**
         * @}
         */

        /**
         * @name Reference quantities
         * @{
         */
        /**
         * Return a reference value typical of the viscosities that appear in
         * this model. This value is not actually used in the material
         * description itself, but is used in scaling variables to the same
         * numerical order of magnitude when solving linear systems.
         * Specifically, the reference viscosity appears in the factor scaling
         * the pressure against the velocity. It is also used in computing
         * dimension-less quantities. You may want to take a look at the
         * Kronbichler, Heister, Bangerth 2012 paper that describes the
         * design of ASPECT for a description of this pressure scaling.
         *
         * @note The reference viscosity should take into account the complete
         * constitutive relationship, defined as the scalar viscosity times the
         * constitutive tensor. In most cases, the constitutive tensor will simply
         * be the identity tensor (this is the default case), but this may become
         * important for material models with anisotropic viscosities, if the
         * constitutive tensor is not normalized.
         */
        virtual double reference_viscosity () const = 0;
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual
        void evaluate (const MaterialModel::MaterialModelInputs<dim> &in,
                       MaterialModel::MaterialModelOutputs<dim> &out) const = 0;
        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);
        /**
         * @}
         */

        /**
         * If this material model can produce additional named outputs
         * that are derived from NamedAdditionalOutputs, create them in here.
         * By default, this does nothing.
          */
        virtual
        void
        create_additional_named_outputs (MaterialModelOutputs &outputs) const;


        /**
         * Fill the additional material model inputs that have been attached
         * by the individual heating or material models in the
         * create_additional_material_model_inputs function.
         * This is done by looping over all material model inputs that have
         * been created and calling their respective member functions.
         */
        virtual
        void
        fill_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &input,
                                              const LinearAlgebra::BlockVector        &solution,
                                              const FEValuesBase<dim>                 &fe_values,
                                              const Introspection<dim>                &introspection) const;

      protected:
        /**
         * A structure that describes how each of the model's
         * output variables (such as viscosity, density, etc) depend
         * on the input variables pressure, temperature, strain rate,
         * and compositional fields.
         *
         * The constructor of this class calls the default
         * constructor of this member variable which in turn
         * initializes the object to invalid values. Derived classes
         * then need to fill it either in their constructor (if they
         * already know the correct dependencies at that time) or
         * at the end of their parse_parameter() functions where
         * they know the correct material parameters they will
         * use.
         */
        NonlinearDependence::ModelDependence model_dependence;
    };

    /**
     * Register a material model so that it can be selected from the parameter
     * file.
     *
     * @param name A string that identifies the material model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this material model wants to read
     * from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this material model.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    void
    register_material_model (const std::string &name,
                             const std::string &description,
                             void (*declare_parameters_function) (ParameterHandler &),
                             Interface<dim> *(*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The material model object returned is not yet initialized and has not
     * read its runtime parameters yet.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    Interface<dim> *
    create_material_model (const std::string &model_name);


    /**
     * A function that reads the name of a model from the parameter object
     * and then returns a pointer to an
     * object that describes this model. Ownership of the pointer is transferred to
     * the caller.
     *
     * The material model object returned is not yet initialized and has not
     * read its runtime parameters yet.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    Interface<dim> *
    create_material_model (ParameterHandler &prm);


    /**
     * Return a string that consists of the names of material models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();


    /**
     * Declare the runtime parameters of the registered material models.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);



    /**
     * For the current plugin subsystem, write a connection graph of all of the
     * plugins we know about, in the format that the
     * programs dot and neato understand. This allows for a visualization of
     * how all of the plugins that ASPECT knows about are interconnected, and
     * connect to other parts of the ASPECT code.
     *
     * @param output_stream The stream to write the output to.
     */
    template <int dim>
    void
    write_plugin_graph (std::ostream &output_stream);



// --------------------- template function definitions ----------------------------------

    template <int dim>
    template <class AdditionalInputType>
    AdditionalInputType *MaterialModelInputs<dim>::get_additional_input()
    {
      for (unsigned int i=0; i<additional_inputs.size(); ++i)
        {
          AdditionalInputType *result = dynamic_cast<AdditionalInputType *> (additional_inputs[i].get());
          if (result)
            return result;
        }
      return nullptr;
    }


    template <int dim>
    template <class AdditionalInputType>
    const AdditionalInputType *MaterialModelInputs<dim>::get_additional_input() const
    {
      for (unsigned int i=0; i<additional_inputs.size(); ++i)
        {
          const AdditionalInputType *result = dynamic_cast<const AdditionalInputType *> (additional_inputs[i].get());
          if (result)
            return result;
        }
      return nullptr;
    }


    template <int dim>
    template <class AdditionalOutputType>
    AdditionalOutputType *MaterialModelOutputs<dim>::get_additional_output()
    {
      for (unsigned int i=0; i<additional_outputs.size(); ++i)
        {
          AdditionalOutputType *result = dynamic_cast<AdditionalOutputType *> (additional_outputs[i].get());
          if (result)
            return result;
        }
      return nullptr;
    }


    template <int dim>
    template <class AdditionalOutputType>
    const AdditionalOutputType *MaterialModelOutputs<dim>::get_additional_output() const
    {
      for (unsigned int i=0; i<additional_outputs.size(); ++i)
        {
          const AdditionalOutputType *result = dynamic_cast<const AdditionalOutputType *> (additional_outputs[i].get());
          if (result)
            return result;
        }
      return nullptr;
    }


    template <int dim>
    void MaterialModelOutputs<dim>::move_additional_outputs_from(MaterialModelOutputs<dim> &other)
    {
      Assert(this->additional_outputs.empty(), ExcMessage("Destination of move needs to be empty!"));
      this->additional_outputs = std::move(other.additional_outputs);
    }


    /**
     * Given a class name, a name, and a description for the parameter file
     * for a material model, register it with the functions that can declare
     * their parameters and create these objects.
     *
     * @ingroup MaterialModels
     */
#define ASPECT_REGISTER_MATERIAL_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_MATERIAL_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::MaterialModel::register_material_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MaterialModel::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::MaterialModel::register_material_model<3>, \
                                name, description); \
  }
  }
}


#endif
