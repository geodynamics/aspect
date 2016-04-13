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


#ifndef __aspect__parameters_h
#define __aspect__parameters_h

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
        IMPES,
        iterated_IMPES,
        iterated_Stokes,
        Stokes_only,
        Advection_only
      };
    };

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

    double                         nonlinear_tolerance;
    bool                           resume_computation;
    double                         start_time;
    double                         CFL_number;
    double                         maximum_time_step;
    bool                           use_artificial_viscosity_smoothing;
    bool                           use_conduction_timestep;
    bool                           convert_to_years;
    std::string                    output_directory;
    double                         surface_pressure;
    double                         adiabatic_surface_temperature;
    unsigned int                   timing_output_frequency;
    bool                           use_direct_stokes_solver;
    double                         linear_stokes_solver_tolerance;
    double                         linear_solver_A_block_tolerance;
    double                         linear_solver_S_block_tolerance;
    unsigned int                   max_nonlinear_iterations;
    unsigned int                   max_nonlinear_iterations_in_prerefinement;
    unsigned int                   n_cheap_stokes_solver_steps;
    double                         temperature_solver_tolerance;
    double                         composition_solver_tolerance;

    /**
     * @}
     */

    /**
     * @name Parameters that have to do with terms in the model
     * @{
     */
    bool                           include_melt_transport;
    double                         melt_transport_threshold;
    double                         radiogenic_heating_rate;
    std::set<types::boundary_id> fixed_temperature_boundary_indicators;
    std::set<types::boundary_id> fixed_composition_boundary_indicators;
    std::set<types::boundary_id> zero_velocity_boundary_indicators;
    std::set<types::boundary_id> tangential_velocity_boundary_indicators;

    /**
     * Map from boundary id to a pair "components", "velocity boundary type",
     * where components is of the format "[x][y][z]" and the velocity type is
     * mapped to one of the plugins of velocity boundary conditions (e.g.
     * "function")
     */
    std::map<types::boundary_id, std::pair<std::string,std::string> > prescribed_velocity_boundary_indicators;

    /**
     * Map from boundary id to a pair "components", "traction boundary type",
     * where components is of the format "[x][y][z]" and the traction type is
     * mapped to one of the plugins of traction boundary conditions (e.g.
     * "function")
     */
    std::map<types::boundary_id, std::pair<std::string,std::string> > prescribed_traction_boundary_indicators;

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
    unsigned int                   min_grid_level;
    std::vector<double>            additional_refinement_times;
    unsigned int                   adaptive_refinement_interval;
    bool                           run_postprocessors_on_initial_refinement;
    /**
     * @}
     */

    /**
     * @name Parameters that have to do with the stabilization of transport
     * equations
     * @{
     */
    unsigned int                   stabilization_alpha;
    double                         stabilization_c_R;
    double                         stabilization_beta;
    double                         discontinuous_penalty;
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
     * @name Parameters that have to do with compositional fields
     * @{
     */
    unsigned int                   n_compositional_fields;
    std::vector<std::string>       names_of_compositional_fields;
    std::vector<unsigned int>      normalized_fields;
    /**
     * @}
     */
    /**
     * @name Parameters that have to do with free surface
     * @{
     */
    bool                           free_surface_enabled;
    std::set<types::boundary_id> free_surface_boundary_indicators;
    /**
     * @}
     */

  };

}
#endif
