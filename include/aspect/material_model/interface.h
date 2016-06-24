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


#ifndef __aspect__material_model_interface_h
#define __aspect__material_model_interface_h

#include <aspect/plugins.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>

namespace aspect
{
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
     * An namespace whose enum members are used in querying the nonlinear
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
        return Dependence((int)d1 | (int)d2);
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
       * Optional reference to the cell that contains these quadrature
       * points. This allows for evaluating properties at the cell vertices
       * and interpolating to the quadrature points, or to query the cell for
       * material ids, neighbors, or other information that is not available
       * solely from the locations. Note that not all calling functions can set
       * this reference. In these cases it will be a NULL pointer, so make sure
       * that your material model either fails with a proper error message
       * or provide an alternative calculation for these cases.
       */
      const typename DoFHandler<dim>::active_cell_iterator *cell;
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
       * Viscosity $\eta$ values at the given positions.
       */
      std::vector<double> viscosities;

      /**
       * Stress-strain "director" tensors at the given positions. This
       * variable can be used to implement exotic rheologies such as
       * anisotropic viscosity.
       *
       * @note The strain rate term in equation (1) of the manual will be
       * multiplied by this tensor *and* the viscosity scalar ($\eta$), as
       * described in the manual secion titled "Constitutive laws". This
       * variable is assigned the rank-four identity tensor by default.
       * This leaves the isotropic constitutive law unchanged if the material
       * model does not explicitly assign a value.
       */
      std::vector<SymmetricTensor<4,dim> > stress_strain_directors;

      /**
       * Density values at the given positions.
       */
      std::vector<double> densities;

      /**
       * Thermal expansion coefficients at the given positions.
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
       * Compressibility at the given positions. The compressibility is given
       * as $\frac 1\rho \frac{\partial\rho}{\partial p}$.
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
      std::vector<std_cxx11::shared_ptr<AdditionalMaterialOutputs<dim> > > additional_outputs;

      /**
       * Given an additional material model output class as explicitly specified
       * template argument, returns a pointer to this additional material model
       * output object if it used in the current simulation.
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
       */
      enum AveragingOperation
      {
        none,
        arithmetic_average,
        harmonic_average,
        geometric_average,
        pick_largest,
        project_to_Q1,
        log_average
      };


      /**
       * Return a string that represents the various averaging options laid
       * out above and that can be used in the claration of an input
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
     * A base class for additional output fields to be added to the
     * MaterialModel::MaterialModelOutputs structure and filled in the
     * MaterialModel::Interface::evaluate() function. The format of the
     * additional quantities defined in derived classes should be the
     * same as for MaterialModel::MaterialModelOutputs.
     */
    template<int dim>
    class AdditionalMaterialOutputs
    {
      public:
        virtual ~AdditionalMaterialOutputs()
        {}

        virtual void average (const MaterialAveraging::AveragingOperation /*operation*/,
                              const FullMatrix<double>  &/*projection_matrix*/,
                              const FullMatrix<double>  &/*expansion_matrix*/)
        {}
    };



    /**
     * A base class for parameterizations of material models. Classes derived
     * from this class will need to implement functions that provide material
     * parameters such as the viscosity, density, etc, typically as a function
     * of position, temperature and pressure at that location.
     *
     * There is two ways to implement a material model and they can not be
     * mixed: Option one is to override all the virtual functions like
     * viscosity(), density(), etc. but not change evaluate().
     *
     * Option two only requires you to override evaluate() and fill the output
     * argument struct instead of implementing the functions viscosity(),
     * density(), etc.. In this case, all other functions are being ignored.
     *
     * The second option is more efficient in general, but it is okay to use
     * option one for simple material models.
     *
     * In all cases, model_dependence values, is_compressible(), reference_viscosity(),
     * reference_density() need to be implemented.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * A typedef to import the MatertialModelInputs name into the current
         * class. This typedef primarily exists as a backward compatibility
         * measure given that the referenced structure used to be a member of
         * the current class.
         */
        typedef MaterialModel::MaterialModelInputs<dim> MaterialModelInputs;
        /**
         * A typedef to import the MatertialModelOutputs name into the current
         * class. This typedef primarily exists as a backward compatibility
         * measure given that the referenced structure used to be a member of
         * the current class.
         */
        typedef MaterialModel::MaterialModelOutputs<dim> MaterialModelOutputs;

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
         * @name Physical parameters used in the basic equations
         * @{
         */
        /**
         * Return the viscosity ratio between disclocation creep and diffusion
         * creep in the case of composite rheology
         */
        virtual double viscosity_ratio (const double      temperature,
                                        const double      pressure,
                                        const std::vector<double>    &compositional_fields,
                                        const SymmetricTensor<2,dim> &strainrate,
                                        const Point<dim> &position) const;


        /**
         * @}
         */

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
         * dimension-less quantities.
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
         * Return the reference density $\rho$. Like the value returned by
         * reference_viscosity(), this value is not actually used in
         * computations but only in postprocessing such as when computing
         * dimension-less quantities.
         */
        virtual double reference_density () const = 0;

        /**
         * Return a reference value for the thermal expansion coefficient
         * $\alpha$. See the thermal_expansion_coefficient() function for a
         * definition of $\alpha$.
         */
        virtual double reference_thermal_expansion_coefficient () const;
        /**
         * @}
         */

        /**
         * @name Auxiliary material properties used for postprocessing
         * @{
         */
        /**
         * Return the p-wave seismic velocity Vp of the model as a function of
         * temperature and pressure.
         *
         * This function is only called in postprocessing. Derived classes do
         * not need to implement it if no useful information is known to
         * compute this quantity, in which case graphical output will simply
         * show an uninformative field of constant value. By default this
         * function returns -1 to indicate that no useful value is
         * implemented.
         */
        virtual
        double
        seismic_Vp (const double      temperature,
                    const double      pressure,
                    const std::vector<double> &compositional_fields,
                    const Point<dim> &position) const;

        /**
         * Return the s-wave seismic velocity Vs of the model as a function of
         * temperature and pressure.
         *
         * This function is only called in postprocessing. Derived classes do
         * not need to implement it if no useful information is known to
         * compute this quantity, in which case graphical output will simply
         * show an uninformative field of constant value. By default this
         * function returns -1 to indicate that no useful value is
         * implemented.
         */
        virtual
        double
        seismic_Vs (const double      temperature,
                    const double      pressure,
                    const std::vector<double> &compositional_fields,
                    const Point<dim> &position) const;
        /**
         * Return the Phase number of the model as a function of temperature
         * and pressure.
         *
         * This function is only called in postprocessing. Derived classes do
         * not need to implement it if no useful information is known to
         * compute this quantity, in which case graphical output will simply
         * show an uninformative field of constant value. By default this
         * function returns 0 to indicate everything is part of the same phase
         */
        virtual
        unsigned int
        thermodynamic_phase (const double      temperature,
                             const double      pressure,
                             const std::vector<double> &compositional_fields) const;
        /**
         * @}
         */

        /**
         * Function to compute the material properties in @p out given the
         * inputs in @p in. If MaterialModelInputs.strain_rate has the length
         * 0, then the viscosity does not need to be computed.
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
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
         * already know the correct dependences at that time) or
         * at the end of their parse_parameter() functions where
         * they know the correct material parameters they will
         * use.
         */
        NonlinearDependence::ModelDependence model_dependence;
    };


    /**
     * This class allows material models written in the past to be used
     * without adapting them to the new interface that requires implementing a
     * function evaluate() for the physical properties. Derive from this
     * helper class instead of Interface and implement the virtual functions
     * viscosity(), etc..
     *
     * Note: do not use this class for new material models, but derive from
     * Interface instead.
     */
    template <int dim>
    class InterfaceCompatibility: public Interface<dim>
    {
      public:
        /**
         * Return the viscosity $\eta$ of the model as a function of
         * temperature, pressure, composition, strain rate, and position.
         *
         * @note The strain rate given as the third argument of this function
         * is computed as $\varepsilon(\mathbf u)=\frac 12 (\nabla \mathbf u +
         * \nabla \mathbf u^T)$, regardless of whether the model is
         * compressible or not. This is relevant since in some other contexts,
         * the strain rate in the compressible case is computed as
         * $\varepsilon(\mathbf u)=\frac 12 (\nabla \mathbf u + \nabla \mathbf
         * u^T) - \frac 13 \nabla \cdot \mathbf u \mathbf 1$.
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const=0;


        /**
         * Return the density $\rho$ of the model as a function of
         * temperature, pressure and position.
         */
        virtual double density (const double      temperature,
                                const double      pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const=0;

        /**
         * Return the compressibility coefficient $\frac 1\rho
         * \frac{\partial\rho}{\partial p}$ of the model as a function of
         * temperature, pressure and position.
         *
         * The compressibility can equivalently be computed as $-\frac 1V
         * \frac{\partial V}{\partial p}$. Note the difference in sign.
         */
        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const=0;

        /**
         * Return the specific heat $C_p$ of the model as a function of
         * temperature, pressure and position.
         */
        virtual double specific_heat (const double      temperature,
                                      const double      pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const=0;

        /**
         * Return the thermal expansion coefficient $\alpha$ of the model,
         * possibly as a function of temperature, pressure and position. The
         * thermal expansion coefficient is defined as $\alpha=-\frac{1}{\rho}
         * \frac{d\rho}{dT}$. Since the density <i>decreases</i> with
         * temperature for almost all models, $\alpha$ is usually positive.
         *
         * The thermal expansion coefficient can equivalently be computed as
         * $\frac 1V \frac{\partial V}{\partial T}$. Note the difference in
         * sign.
         */
        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const=0;

        /**
         * Return the product of the change in entropy across phase
         * transitions, the pressure derivative of the phase function (if this
         * is the pressure derivative) or the product of the former two and
         * the Clapeyron slope (if this is the temperature derivative). The
         * entropy change across a phase transition can be calculated as
         * $\frac{\gamma \Delta\rho}{\rho_\text{light} \rho_\text{heavy}}$.
         * $\gamma$ is the Clapeyron slope of the phase transition,
         * $\Delta\rho$ is the density jump across the phase transition,
         * $\rho_\text{light}$ is the density of the light material (above the
         * phase transition) and $\rho_\text{heavy}$ the density of the heavy
         * material (below the phase transition). The phase function hat
         * values ranging from 0 to 1 indicating which percentage of the
         * material has already undergone the phase transition. Its argument
         * is usually the excess pressure $\pi = p - p_0 - \gamma T$, where
         * $p_0$ is the zero-degree transition pressure.
         *
         * This function has a default implementation that sets the entropy
         * gradient to zero (assuming no phase changes).
         */
        virtual double entropy_derivative (const double      temperature,
                                           const double      pressure,
                                           const std::vector<double> &compositional_fields,
                                           const Point<dim> &position,
                                           const NonlinearDependence::Dependence dependence) const;

        /**
         * Return the change in the compositional field compositional_variable
         * due to reactions between different compositional fields. It is
         * assumed that there is always an equilibrium between the
         * compositional fields (because the time scale of reactions is
         * normally much shorter than that of convection), so the quantity
         * returned by this function is an actual change in the amount of
         * material, which is added to or subtracted from the current value of
         * the compositional field, and NOT a reaction rate. The idea is, that
         * in dependence of temperature, pressure, position and the
         * compositional fields themselves an equilibrium can be calculated,
         * and the difference between the current value and the equilibrium
         * can be added to the respective compositional field.
         *
         * For mass conservation it should ALWAYS be checked that what is
         * subtracted from one field is added to another field (and the other
         * way round) and that one never subtracts more than the actual value
         * of a field (so it does not get negative).
         *
         * This function has a default implementation that sets the reaction
         * term to zero (assuming no reactions).
         *
         * @note In cases where one has slow chemical reactions (or cases
         * where compositional fields are used to track quantities different
         * than actual compositions, for example accumulated strains in damage
         * models), models are formulated as differential equations with right
         * hand sides, not as instantaneous equations. In such cases, the
         * reaction terms (i.e., the incremental additions to the previous
         * state) are usually of the form reaction rate times time step size.
         * To implement something like this, derive your material model from
         * SimulatorAccess so you can query the time step used by the
         * simulator in order to compute the reaction increment.
         */
        virtual double reaction_term (const double      temperature,
                                      const double      pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position,
                                      const unsigned int compositional_variable) const;

        /**
         * Return the thermal conductivity $k$ of the model as a function of
         * temperature, pressure and position. The units of $k$ are
         * $\textrm{W} / \textrm{m} / \textrm{K}$ in 3d, and $\textrm{W} /
         * \textrm{K}$ in 2d. This is easily see by considering that $k$ is
         * the heat flux density (i.e., Watts per unit area perpendicular to
         * the heat flux direction) per unit temperature gradient (i.e.,
         * Kelvin per meter). The unit area has units $m^2$ in 3d, but only
         * $m$ in 2d, yielding the stated units for $k$.
         *
         * Note that the thermal <i>conductivity</i> $k$ is related to the
         * thermal <i>diffusivity</i> $\kappa$ as $k = \kappa \rho c_p$. In
         * essence, the conductivity relates to the question of how thermal
         * energy diffuses whereas the diffusivity relates to the question of
         * how the temperature diffuses. $\kappa$ has units
         * $\textrm{m}^2/\textrm{s}$.
         */
        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const=0;


        /**
         * The evaluate() function is implemented to call the individual
         * functions in this class, so there is no need to implement this in
         * your material model derived from InterfaceCompatibility.
         * @param in
         * @param out
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;
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
      return NULL;
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
      return NULL;
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
