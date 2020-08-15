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


#ifndef _aspect_parameters_h
#define _aspect_parameters_h

#include <deal.II/base/parameter_handler.h>

#include <aspect/global.h>
#include <aspect/material_model/interface.h>


namespace aspect
{
  using namespace dealii;

  // forward declaration:
  namespace GeometryModel
  {
    template <int dim>
    class Interface;
  }

  /**
   * A structure that holds run-time parameters. These parameters are all
   * declared for the ParameterHandler class in the declare_parameters()
   * member function, and read in the parse_parameters() function.
   *
   * Each of the member variables of this class corresponds to a parameter
   * declared for the ParameterHandler class. Rather than duplicating the
   * documentation of each of these parameters for the member variables here,
   * please refer to the documentation of run-time parameters in the ASPECT
   * manual for more information.
   *
   * @ingroup Simulator
   */
  template <int dim>
  struct Parameters
  {

    /**
     * A structure that contains enum values that identify the nonlinear
     * solver in use.
     */
    struct NonlinearSolver
    {
      enum Kind
      {
        single_Advection_single_Stokes,
        iterated_Advection_and_Stokes,
        single_Advection_iterated_Stokes,
        no_Advection_iterated_Stokes,
        no_Advection_single_Stokes,
        no_Advection_iterated_defect_correction_Stokes,
        single_Advection_iterated_defect_correction_Stokes,
        iterated_Advection_and_defect_correction_Stokes,
        iterated_Advection_and_Newton_Stokes,
        single_Advection_iterated_Newton_Stokes,
        single_Advection_no_Stokes,
        first_timestep_only_single_Stokes,
        no_Advection_no_Stokes
      };
    };

    static
    bool
    is_defect_correction(const typename NonlinearSolver::Kind &input)
    {
      return input == NonlinearSolver::iterated_Advection_and_Newton_Stokes ||
             input == NonlinearSolver::single_Advection_iterated_Newton_Stokes ||
             input == NonlinearSolver::no_Advection_iterated_defect_correction_Stokes ||
             input == NonlinearSolver::single_Advection_iterated_defect_correction_Stokes ||
             input == NonlinearSolver::iterated_Advection_and_defect_correction_Stokes
             ?
             true
             :
             false;
    }

    /**
     * @brief The NullspaceRemoval struct
     */
    struct NullspaceRemoval
    {
      enum Kind
      {
        none = 0,
        net_translation_x = 0x1,
        net_translation_y = 0x2,
        net_translation_z = 0x4,
        net_translation   = 0x1+0x2+0x4,
        linear_momentum_x = 0x8,
        linear_momentum_y = 0x10,
        linear_momentum_z = 0x20,
        linear_momentum   = 0x8+0x10+0x20,
        net_rotation      = 0x40,
        angular_momentum  = 0x80
      };
    };

    /**
     * A struct that describes the available methods to solve
     * advected fields. This type is at the moment only used to determine how
     * to advect each compositional field -- or what else to do with it if
     * it doesn't satisfy an advection equation, for example if the
     * compositional field just contains data interpolated from
     * particles in each time step (`particles`) or if it contains
     * data interpolated from other sources such as material outputs
     * (`prescribed_field`), neither of which is further advected.
     */
    struct AdvectionFieldMethod
    {
      enum Kind
      {
        fem_field,
        particles,
        volume_of_fluid,
        static_field,
        fem_melt_field,
        prescribed_field,
        prescribed_field_with_diffusion
      };
    };

    /**
     * A struct that contains information about which
     * formulation of the basic equations should be solved,
     * i.e. which terms to consider, which terms to neglect, and
     * which to simplify in different ways.
     */
    struct Formulation
    {
      /**
       * This enum lists available formulations that
       * determine the combined approximations made in
       * all solved equations. 'Custom' allows to set
       * approximations individually for single equations.
       */
      enum Kind
      {
        boussinesq_approximation,
        anelastic_liquid_approximation,
        isentropic_compression,
        custom
      };

      /**
       * This function translates an input string into the
       * available enum options.
       */
      static
      Kind
      parse(const std::string &input)
      {
        if (input == "isentropic compression")
          return Formulation::isentropic_compression;
        else if (input == "anelastic liquid approximation")
          return Formulation::anelastic_liquid_approximation;
        else if (input == "Boussinesq approximation")
          return Formulation::boussinesq_approximation;
        else if (input == "custom")
          return Formulation::custom;
        else
          AssertThrow(false, ExcNotImplemented());

        return Formulation::Kind();
      }

      /**
       * This struct contains information about the approximation
       * made to the term containing the gradient of the density
       * in the mass conservation equation. The different possible
       * ways to approximate this term are described in the manual.
       */
      struct MassConservation
      {
        /**
         * This enum lists available approximations to the
         * density gradient term of the mass conservation equation.
         */
        enum Kind
        {
          isentropic_compression,
          hydrostatic_compression,
          reference_density_profile,
          implicit_reference_density_profile,
          incompressible,
          projected_density_field,
          ask_material_model
        };

        /**
         * This function translates an input string into the
         * available enum options.
         */
        static Kind
        parse(const std::string &input)
        {
          if (input == "isentropic compression")
            return Formulation::MassConservation::isentropic_compression;
          else if (input == "hydrostatic compression")
            return Formulation::MassConservation::hydrostatic_compression;
          else if (input == "reference density profile")
            return Formulation::MassConservation::reference_density_profile;
          else if (input == "implicit reference density profile")
            return Formulation::MassConservation::implicit_reference_density_profile;
          else if (input == "incompressible")
            return Formulation::MassConservation::incompressible;
          else if (input == "projected density field")
            return Formulation::MassConservation::projected_density_field;
          else if (input == "ask material model")
            return Formulation::MassConservation::ask_material_model;
          else
            AssertThrow(false, ExcNotImplemented());

          return Formulation::MassConservation::Kind();
        }
      };

      /**
       * This struct contains information about the approximation
       * made to temperature equation. The different possible
       * ways to approximate this term are described in the manual.
       */
      struct TemperatureEquation
      {
        /**
         * This enum lists available approximations to the
         * density in the temperature equation.
         */
        enum Kind
        {
          real_density,
          reference_density_profile
        };

        /**
         * This function translates an input string into the
         * available enum options.
         */
        static
        Kind
        parse(const std::string &input)
        {
          if (input == "real density")
            return Formulation::TemperatureEquation::real_density;
          else if (input == "reference density profile")
            return Formulation::TemperatureEquation::reference_density_profile;
          else
            AssertThrow(false, ExcNotImplemented());

          return Formulation::TemperatureEquation::Kind();
        }
      };
    };

    /**
     * A struct to provide the setting for "Stabilization method"
     */
    struct AdvectionStabilizationMethod
    {
      enum Kind
      {
        entropy_viscosity,
        supg
      };

      /**
       * This function translates an input string into the
       * available enum options.
       */
      static
      Kind
      parse(const std::string &input)
      {
        if (input == "entropy viscosity")
          return entropy_viscosity;
        else if (input == "SUPG")
          return supg;
        else
          AssertThrow(false, ExcNotImplemented());

        return Kind();
      }

      static std::string get_options_string()
      {
        return "entropy viscosity|SUPG";
      }
    };

    /**
     * This enum represents the different choices for the linear solver
     * for the Stoke system. See @p stokes_solver_type.
     */
    struct StokesSolverType
    {
      enum Kind
      {
        block_amg,
        direct_solver,
        block_gmg
      };

      static const std::string pattern()
      {
        return "block AMG|direct solver|block GMG";
      }

      static Kind
      parse(const std::string &input)
      {
        if (input == "block AMG")
          return block_amg;
        else if (input == "direct solver")
          return direct_solver;
        else if (input == "block GMG")
          return block_gmg;
        else
          AssertThrow(false, ExcNotImplemented());

        return Kind();
      }
    };

    /**
     * This enum represents the different choices for the Krylov method
     * used in the cheap GMG Stokes solve.
     */
    struct StokesKrylovType
    {
      enum Kind
      {
        gmres,
        idr_s
      };

      static const std::string pattern()
      {
        return "GMRES|IDR(s)";
      }

      static Kind
      parse(const std::string &input)
      {
        if (input == "GMRES")
          return gmres;
        else if (input == "IDR(s)")
          return idr_s;
        else
          AssertThrow(false, ExcNotImplemented());

        return Kind();
      }
    };

    /**
     * Constructor. Fills the values of member functions from the given
     * parameter object.
     *
     * @param prm The parameter object that has previously been filled with
     * content by reading an input file.
     *
     * @param mpi_communicator The MPI communicator we will use for this
     * simulation. We need this when calling parse_parameters() so that we can
     * verify some of the input arguments.
     */
    Parameters (ParameterHandler &prm,
                MPI_Comm mpi_communicator);

    /**
     * Declare the run-time parameters this class takes, and call the
     * respective <code>declare_parameters</code> functions of the namespaces
     * that describe geometries, material models, etc.
     *
     * @param prm The object in which the run-time parameters are to be
     * declared.
     */
    static
    void declare_parameters (ParameterHandler &prm);

    /**
     * Read run-time parameters from an object that has previously parsed an
     * input file. This reads all parameters that do not require knowledge of
     * the geometry model we use. There is a separate function
     * parse_geometry_dependent_parameters() that is called as soon as the
     * geometry object has been created and that can translate between the
     * symbolic names for boundary components that the geometry model
     * publishes and the boundary indicators used internally.
     *
     * @param prm The object from which to obtain the run-time parameters.
     *
     * @param mpi_communicator The MPI communicator we will use for this
     * simulation. We need this when calling parse_parameters() so that we can
     * verify some of the input arguments.
     */
    void parse_parameters (ParameterHandler &prm,
                           const MPI_Comm mpi_communicator);

    /**
     * Read those run-time parameters from a ParameterHandler object that
     * depend on knowing which geometry object we use. This function
     * complements parse_parameters() but is only called once the geometry
     * object has been created. This function is separate because we allow the
     * use of symbolic names in defining which boundary components have which
     * boundary conditions, and the names one can specify there are not
     * available until after the geometry object has been created.
     *
     * This function is called from the GeometryModel::create_geometry()
     * function.
     *
     * @param prm The object from which to obtain the run-time parameters.
     * @param geometry_model The geometry model that provides boundary names
     * etc.
     */
    void parse_geometry_dependent_parameters (ParameterHandler &prm,
                                              const GeometryModel::Interface<dim> &geometry_model);

    /**
     * @name Global parameters
     * @{
     */
    typename NonlinearSolver::Kind nonlinear_solver;

    typename AdvectionStabilizationMethod::Kind advection_stabilization_method;
    double                         nonlinear_tolerance;
    bool                           resume_computation;
    double                         start_time;
    double                         CFL_number;
    double                         maximum_time_step;
    double                         maximum_relative_increase_time_step;
    double                         maximum_first_time_step;
    bool                           use_artificial_viscosity_smoothing;
    bool                           use_conduction_timestep;
    bool                           convert_to_years;
    std::string                    output_directory;
    double                         surface_pressure;
    double                         adiabatic_surface_temperature;
    unsigned int                   timing_output_frequency;
    unsigned int                   max_nonlinear_iterations;
    unsigned int                   max_nonlinear_iterations_in_prerefinement;
    bool                           use_operator_splitting;
    std::string                    world_builder_file;

    /**
     * @}
     */

    /**
     * @name section: Solver parameters
     * @{
     */
    double                         temperature_solver_tolerance;
    double                         composition_solver_tolerance;

    // subsection: Advection solver parameters
    unsigned int                   advection_gmres_restart_length;

    // subsection: Stokes solver parameters
    bool                           use_direct_stokes_solver;
    typename StokesSolverType::Kind stokes_solver_type;
    typename StokesKrylovType::Kind stokes_krylov_type;
    unsigned int                    idr_s_parameter;

    double                         linear_stokes_solver_tolerance;
    unsigned int                   n_cheap_stokes_solver_steps;
    unsigned int                   n_expensive_stokes_solver_steps;
    double                         linear_solver_A_block_tolerance;
    bool                           use_full_A_block_preconditioner;
    double                         linear_solver_S_block_tolerance;
    unsigned int                   stokes_gmres_restart_length;

    // subsection: AMG parameters
    std::string                    AMG_smoother_type;
    unsigned int                   AMG_smoother_sweeps;
    double                         AMG_aggregation_threshold;
    bool                           AMG_output_details;

    // subsection: Operator splitting parameters
    double                         reaction_time_step;
    unsigned int                   reaction_steps_per_advection_step;

    // subsection: Diffusion solver parameters
    double                         diffusion_length_scale;

    /**
     * @}
     */

    /**
     * @name Formulation settings
     * @{
     */

    /**
     * This variable determines which of the several ways to formulate the
     * equations ASPECT will solve.
     * Common formulations are the Boussinesq or Anelastic Liquid
     * Approximations (BA, ALA). ASPECT's original formulation is termed
     * 'isentropic compression'. 'Custom' allows
     * to set the approximations individually per equation.
     */
    typename Formulation::Kind formulation;

    /**
     * Determines how to formulate the mass conservation equation in ASPECT.
     * Common approximations are 'incompressible' or 'reference density profile'.
     * ASPECT's original formulation is termed 'isentropic compression'. See the
     * manual for more details about the individual terms.
     */
    typename Formulation::MassConservation::Kind formulation_mass_conservation;

    /**
     * Determines how to formulate the density in the temperature equation
     * in ASPECT. Possible approximations are 'reference density profile' or
     * 'real density'.
     */
    typename Formulation::TemperatureEquation::Kind formulation_temperature_equation;

    /**
     * This variable determines whether additional terms related to elastic forces
     * are added to the Stokes equation.
     */
    bool                           enable_elasticity;

    /**
     * @}
     */

    /**
     * @name Parameters that have to do with terms in the model
     * @{
     */
    bool                           include_melt_transport;
    bool                           enable_additional_stokes_rhs;
    bool                           enable_prescribed_dilation;

    /**
     * Map from boundary id to a pair "components", "traction boundary type",
     * where components is of the format "[x][y][z]" and the traction type is
     * mapped to one of the plugins of traction boundary conditions (e.g.
     * "function")
     */
    std::map<types::boundary_id, std::pair<std::string,std::string> > prescribed_traction_boundary_indicators;

    /**
     * A set of boundary ids on which the boundary_heat_flux objects
     * will be applied.
     */
    std::set<types::boundary_id> fixed_heat_flux_boundary_indicators;

    /**
     * Selection of operations to perform to remove nullspace from velocity
     * field.
     */
    typename NullspaceRemoval::Kind nullspace_removal;
    /**
     * @}
     */

    /**
     * @name Parameters that have to do with mesh refinement
     * @{
     */
    unsigned int                   initial_global_refinement;
    unsigned int                   initial_adaptive_refinement;
    double                         refinement_fraction;
    double                         coarsening_fraction;
    bool                           adapt_by_fraction_of_cells;
    unsigned int                   min_grid_level;
    std::vector<double>            additional_refinement_times;
    unsigned int                   adaptive_refinement_interval;
    bool                           skip_solvers_on_initial_refinement;
    bool                           skip_setup_initial_conditions_on_initial_refinement;
    bool                           run_postprocessors_on_initial_refinement;
    bool                           run_postprocessors_on_nonlinear_iterations;
    /**
     * @}
     */

    /**
     * @name Parameters that have to do with the stabilization of transport
     * equations
     * @{
     */
    unsigned int                   stabilization_alpha;
    std::vector<double>            stabilization_c_R;
    std::vector<double>            stabilization_beta;
    double                         stabilization_gamma;
    double                         discontinuous_penalty;
    bool                           use_limiter_for_discontinuous_temperature_solution;
    bool                           use_limiter_for_discontinuous_composition_solution;
    double                         global_temperature_max_preset;
    double                         global_temperature_min_preset;
    std::vector<double>            global_composition_max_preset;
    std::vector<double>            global_composition_min_preset;
    /**
     * @}
     */

    /**
     * @name Parameters that have to do with checkpointing
     * @{
     */
    int                            checkpoint_time_secs;
    int                            checkpoint_steps;
    /**
     * @}
     */

    /**
     * @name Parameters that have to do with spatial discretizations
     * @{
     */
    unsigned int                   stokes_velocity_degree;
    bool                           use_locally_conservative_discretization;
    bool                           use_equal_order_interpolation_for_stokes;
    bool                           use_discontinuous_temperature_discretization;
    bool                           use_discontinuous_composition_discretization;
    unsigned int                   temperature_degree;
    unsigned int                   composition_degree;
    std::string                    pressure_normalization;
    MaterialModel::MaterialAveraging::AveragingOperation material_averaging;

    /**
     * @}
     */

    /**
     * @name Parameters that have to do with the temperature field
     * @{
     */

    typename AdvectionFieldMethod::Kind temperature_method;

    /**
     * @}
     */

    /**
     * @name Parameters that have to do with compositional fields
     * @{
     */
    unsigned int                   n_compositional_fields;
    std::vector<std::string>       names_of_compositional_fields;

    /**
     * A vector that contains the advection field method for every compositional
     * field. Consequently the vector has n_compositional_fields entries.
     */
    std::vector<typename AdvectionFieldMethod::Kind> compositional_field_methods;

    /**
     * Map from compositional index to a pair "particle property", "component",
     * where particle property is a string that can be mapped to one of the
     * particle property plugins.
     * Component denotes which component of the particle property is to be
     * mapped in case there are several. Therefore, it is optional to specify
     * the component and it is of the format "[0][1][2]". In case no component
     * is specified it defaults to 0.
     */
    std::map<unsigned int, std::pair<std::string,unsigned int> > mapped_particle_properties;

    std::vector<unsigned int>      normalized_fields;
    /**
     * @}
     */
    /**
     * @name Parameters that have to do with mesh deformation
     * @{
     */
    bool                           mesh_deformation_enabled;
    /**
     * @}
     */

    /**
     * @name Parameters that have to do with volume of fluid calculations
     * @{
     */
    bool                           volume_of_fluid_tracking_enabled;
    /**
     * @}
     */


  };

}
#endif
