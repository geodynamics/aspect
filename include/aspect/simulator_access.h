/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__simulator_access_h
#define __aspect__simulator_access_h

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q.h>

#include <aspect/global.h>
#include <aspect/introspection.h>
#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/initial_conditions/interface.h>
#include <aspect/compositional_initial_conditions/interface.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/mesh_refinement/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/heating_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>



namespace aspect
{
  using namespace dealii;

  // forward declaration
  template <int> class Simulator;


  /**
   * SimulatorAccess is base class for different plugins like postprocessors.
   * It provides access to the various variables of the main class that
   * plugins may want to use in their evaluations, such as solution vectors,
   * the current time, time step sizes, material models, or the triangulations
   * and DoFHandlers that correspond to solutions.
   *
   * This class is the interface between plugins and the main simulator class.
   * Using this insulation layer, the plugins need not know anything about the
   * internal details of the simulation class.
   *
   * Every Postprocessor is required to derive from SimulatorAccess. It is
   * optional for other plugins like MaterialModel, GravityModel, etc..
   *
   * Since the functions providing access to details of the simulator class
   * are meant to be used only by derived classes of this class (rather than
   * becoming part of the public interface of these classes), the functions of
   * this class are made @p protected.
   *
   * @ingroup Simulator
   */
  template <int dim>
  class SimulatorAccess
  {
    public:
      /**
       * Destructor. Does nothing but is virtual so that derived classes
       * destructors are also virtual.
       */
      virtual
      ~SimulatorAccess ();

      /**
       * Initialize this class for a given simulator. This function is marked
       * as virtual so that derived classes can do something upon
       * initialization as well, for example look up and cache data; derived
       * classes should call this function from the base class as well,
       * however.
       *
       * @param simulator A reference to the main simulator object.
       */
      virtual void initialize (const Simulator<dim> &simulator);

    protected:
      /** @name Accessing variables that identify overall properties of the simulator */
      /** @{ */

      /**
       * Return a reference to an introspection object that describes overall
       * properties of the simulator. In particular, it provides symbolic
       * names for extractors and component masks for each variable, etc, and
       * thereby reduces the need for implicit knowledge throughout the code
       * base.
       */
      const Introspection<dim> &
      introspection () const;

      /**
       * Returns a reference to the Simulator itself. Note that you can not
       * access any members or functions of the Simulator. This function
       * exists so that any class with SimulatorAccess can create other
       * objects with SimulatorAccess (because initializing them requires a
       * reference to the Simulator).
       */
      const Simulator<dim> &
      get_simulator() const;

      /**
       * Return the MPI communicator for this simulation.
       */
      MPI_Comm
      get_mpi_communicator () const;

      /**
       * Return a reference to the stream object that only outputs something
       * on one processor in a parallel program and simply ignores output put
       * into it on all other processors.
       */
      const ConditionalOStream &
      get_pcout () const;

      /**
       * Return the current simulation time in seconds.
       */
      double get_time () const;

      /**
       * Return the size of the last time step.
       */
      double
      get_timestep () const;

      /**
       * Return the current number of a time step.
       */
      unsigned int
      get_timestep_number () const;

      /**
       * Return a reference to the triangulation in use by the simulator
       * object.
       */
      const parallel::distributed::Triangulation<dim> &
      get_triangulation () const;

      /**
       * Return the global volume of the computational domain.
       */
      double
      get_volume () const;

      /**
       * Return a reference to the mapping used to describe the boundary of
       * the domain.
       */
      const Mapping<dim> &
      get_mapping () const;

      /**
       * Return the directory specified in the input parameter file to be the
       * place where output files are to be placed. The string is terminated
       * by a directory separator (i.e., '/').
       */
      std::string
      get_output_directory () const;

      /**
       * Return whether we use the adiabatic heating term.
       */
      bool
      include_adiabatic_heating () const;

      /**
       * Return whether we use the latent heat term.
       */
      bool
      include_latent_heat () const;

      /**
       * Return the stokes velocity degree.
       */
      int
      get_stokes_velocity_degree () const;


      /**
       * Return the adiabatic surface temperature.
       */
      double
      get_adiabatic_surface_temperature () const;

      /**
       * Return the adiabatic surface pressure.
       */
      double
      get_surface_pressure () const;

      /**
       * Return whether things like velocities should be converted from the
       * seconds in the MKS system to years. The value of this flag is set by
       * the corresponding entry in the input parameter file.
       */
      bool
      convert_output_to_years () const;

      /**
       * Return the number of compositional fields specified in the input
       * parameter file that will be advected along with the flow field.
       */
      unsigned int
      n_compositional_fields () const;

      /**
       * Compute the error indicators in the same way they are normally used
       * for mesh refinement. The mesh is not refined when doing so, but the
       * indicators can be used when generating graphical output to check why
       * mesh refinement is proceeding as it is.
       */
      void
      get_refinement_criteria(Vector<float> &estimated_error_per_cell) const;

      /**
       * Returns the entropy viscosity on each locally owned cell as it is
       * used to stabilize the temperature equation.
       */
      void
      get_artificial_viscosity(Vector<float> &viscosity_per_cell) const;

      /** @} */


      /** @name Accessing variables that identify the solution of the problem */
      /** @{ */


      /**
       * Return a reference to the vector that has the current solution of the
       * entire system, i.e. the velocity and pressure variables as well as
       * the temperature.  This vector is associated with the DoFHandler
       * object returned by get_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_solution () const;

      /**
       * Return a reference to the vector that has the solution of the entire
       * system at the previous time step. This vector is associated with the
       * DoFHandler object returned by get_stokes_dof_handler().
       *
       * @note In general the vector is a distributed vector; however, it
       * contains ghost elements for all locally relevant degrees of freedom.
       */
      const LinearAlgebra::BlockVector &
      get_old_solution () const;

      /**
       * Return a reference to the DoFHandler that is used to discretize the
       * variables at the current time step.
       */
      const DoFHandler<dim> &
      get_dof_handler () const;

      /**
       * Return a reference to the finite element that the DoFHandler that is
       * used to discretize the variables at the current time step is built
       * on. This is the finite element for the entire, couple problem, i.e.,
       * it contains sub-elements for velocity, pressure, temperature and all
       * other variables in this problem (e.g., compositional variables, if
       * used in this simulation).
       */
      const FiniteElement<dim> &
      get_fe () const;

      /**
       * Fill the argument with a set of depth averages of the current
       * temperature field. The function fills a vector that contains average
       * field values over slices of the domain of same depth.
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_temperature(std::vector<double> &values) const;

      /**
       * Fill the argument with a set of depth averages of the current
       * compositional fields. See get_depth_average_temperature.
       *
       * @param composition_index The index of the compositional field whose
       * matrix we want to assemble (0 <= composition_index < number of
       * compositional fields in this problem).
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_composition(const unsigned int composition_index,
                                    std::vector<double> &values) const;

      /**
       * Compute a depth average of the current viscosity
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_viscosity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current velocity magnitude
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_velocity_magnitude(std::vector<double> &values) const;

      /**
       * Compute a depth average of the current sinking velocity
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_sinking_velocity(std::vector<double> &values) const;

      /**
       * Compute a depth average of the seismic shear wave speed: Vs
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_Vs(std::vector<double> &values) const;

      /**
       * Compute a depth average of the seismic pressure wave speed: Vp
       *
       * @param values The output vector of depth averaged values. The
       * function takes the pre-existing size of this vector as the number of
       * depth slices.
       */
      void
      get_depth_average_Vp(std::vector<double> &values) const;
      /** @} */


      /** @name Accessing variables that identify aspects of the simulation */
      /** @{ */

      /**
       * Return a pointer to the material model to access functions like
       * density().
       */
      const MaterialModel::Interface<dim> &
      get_material_model () const;

      /**
       * Return a pointer to the gravity model description.
       */
      const GravityModel::Interface<dim> &
      get_gravity_model () const;

      /**
       * Return a pointer to the geometry model.
       */
      const GeometryModel::Interface<dim> &
      get_geometry_model () const;


      /**
       * Return a pointer to the object that describes the adiabatic
       * conditions.
       */
      const AdiabaticConditions::Interface<dim> &
      get_adiabatic_conditions () const;

      /**
       * Return a pointer to the object that describes the temperature
       * boundary values.
       */
      const BoundaryTemperature::Interface<dim> &
      get_boundary_temperature () const;

      /**
       * Return a pointer to the object that describes the temperature initial
       * values.
       */
      const InitialConditions::Interface<dim> &
      get_initial_conditions () const;


      /**
       * Return a pointer to the object that describes the composition initial
       * values.
       */
      const CompositionalInitialConditions::Interface<dim> &
      get_compositional_initial_conditions () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a fixed temperature.
       */
      const std::set<types::boundary_id> &
      get_fixed_temperature_boundary_indicators () const;

      /**
       * Return a set of boundary indicators that describes which of the
       * boundaries have a fixed composition.
       */
      const std::set<types::boundary_id> &
      get_fixed_composition_boundary_indicators () const;

      /**
       * Return the map of prescribed_velocity_boundary_conditions
       */
      const std::map<types::boundary_id,std_cxx1x::shared_ptr<VelocityBoundaryConditions::Interface<dim> > >
      get_prescribed_velocity_boundary_conditions () const;

      /**
       * Return a pointer to the heating model.
       */
      const HeatingModel::Interface<dim> &
      get_heating_model () const;

      /**
       * A convenience function that copies the values of the compositional
       * fields at the quadrature point q given as input parameter to the
       * output vector composition_values_at_q_point.
       */
      static
      void
      get_composition_values_at_q_point (const std::vector<std::vector<double> > &composition_values,
                                         const unsigned int                      q,
                                         std::vector<double>                    &composition_values_at_q_point);


      /**
       * Find a pointer to a certain postprocessor, if not return a NULL
       * pointer.
       */
      template <typename PostprocessorType>
      PostprocessorType *
      find_postprocessor () const;

      /** @} */

    private:
      /**
       * A pointer to the simulator object to which we want to get access.
       */
      const Simulator<dim> *simulator;
  };

  template <int dim>
  template <typename PostprocessorType>
  inline
  PostprocessorType *
  SimulatorAccess<dim>::find_postprocessor () const
  {
    return simulator->postprocess_manager.template find_postprocessor<PostprocessorType>();
  }
}


#endif
