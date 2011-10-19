//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------
#ifndef __aspect__postprocess_base_h
#define __aspect__postprocess_base_h


#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/dofs/dof_handler.h>


namespace aspect
{
  template <int dim> class Simulator;


  /**
   * A namespace for everything to do with postprocessing solutions every time step
   * or every few time steps.
   *
   * @ingroup Postprocessing
   **/
  namespace Postprocess
  {

    /**
     * This class declares the public interface of postprocessors. Postprocessors
     * must implement a function that can be called at the end of each time step
     * to evaluate the current solution, as well as functions that save the state
     * of the object and restore it (for checkpoint/restart capabilities).
     * 
     * Access to the data of the simulator is granted by the @p protected member functions
     * of the SimulatorAccessor class, i.e., classes implementing this interface will
     * in general want to derive from both this Interface class as well as from the
     * SimulatorAccess class.
     * 
     * @ingroup Postprocessing
     */
    template <int dim>
    class Interface
    {
      public:
	/**
	 * Execute this postprocessor. Derived classes will implement this function
	 * to do whatever they want to do to evaluate the solution at the current
	 * time step.
	 *
	 * @param statistics An object that contains statistics that are collected
	 * throughout the simulation and that will be written to an output file at
	 * the end of each time step. Postprocessors may deposit data in these
	 * tables for later visualization or further processing.
	 * 
	 * @return A pair of strings that will be
	 * printed to the screen after running the postprocessor in two columns;
	 * typically the first column contains a description of what the data is
	 * and the second contains a numerical value of this data. If there is
	 * nothing to print, simply return two empty strings.
	 **/
	virtual 
	std::pair<std::string,std::string> 
	execute (TableHandler &statistics) = 0;
    };

    
    
    /**
     * Base class for postprocessors. This class provides access to
     * the various variables of the main class that postprocessors may want to use
     * in their evaluations, such as solution vectors, the current time, time step
     * sizes, or the triangulations and DoFHandlers that correspond to solutions.
     *
     * This class is the interface between postprocessors and the main simulator
     * class. Using this insulation layer, postprocessors need not know anything
     * about the internal details of the simulation class.
     *
     * Since the functions providing access to details of the simulator class
     * are meant to be used only by derived classes of this class (rather
     * than becoming part of the public interface of these classes), the
     * functions of this class are made @p protected.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class SimulatorAccess
    {
      public:
	/**
	 * Constructor.
	 *
	 * @param simulator A reference to the main simulator object to which the
	 * postprocessor implemented in the derived class should be applied.
	 **/
	SimulatorAccess (const Simulator<dim> &simulator);
	
    protected:
        /** @name Accessing variables that identify overall properties of the simulator */
        /** @{ */

        /**
         * Return the current simulation time.
         */
        double get_time () const;

        /**
         * Return the current number of a time step.
         */
        double get_timestep_number () const;
        /** @} */


        /** @name Accessing variables that identify the solution of the Stokes problem */
        /** @{ */


        /**
         * Return a reference to the vector that has the current solution
         * of the Stokes system, i.e. the velocity and pressure variables.
         * This vector is associated with the DoFHandler object returned by
         * get_stokes_dof_handler().
         *
         * @note In general the vector is a distributed vector; however, it
         * contains ghost elements for all locally relevant degrees of freedom.
         */
        const TrilinosWrappers::Vector &
        get_stokes_solution () const;

        /**
         * Return a reference to the vector that has the solution
         * of the Stokes system at the previous time step.
         * This vector is associated with the DoFHandler object returned by
         * get_stokes_dof_handler().
         *
         * @note In general the vector is a distributed vector; however, it
         * contains ghost elements for all locally relevant degrees of freedom.
         */
        const TrilinosWrappers::Vector &
        get_old_stokes_solution () const;

        /**
         * Return a reference to the DoFHandler that is used to discretize
         * the Stokes system of velocity and pressure.
         */
        const parallel::distributed::DoFhandler<dim> &
        get_stokes_dof_handler () const;
        /** @} */


        /** @name Accessing variables that identify the solution of the temperature problem */
        /** @{ */


        /**
         * Return a reference to the vector that has the current solution
         * of the temperature system.
         * This vector is associated with the DoFHandler object returned by
         * get_temperature_dof_handler().
         *
         * @note In general the vector is a distributed vector; however, it
         * contains ghost elements for all locally relevant degrees of freedom.
         */
        const TrilinosWrappers::Vector &
        get_temperature_solution () const;

        /**
         * Return a reference to the vector that has the solution
         * of the temperature system at the previous time step.
         * This vector is associated with the DoFHandler object returned by
         * get_temperature_dof_handler().
         *
         * @note In general the vector is a distributed vector; however, it
         * contains ghost elements for all locally relevant degrees of freedom.
         */
        const TrilinosWrappers::Vector &
        get_old_temperature_solution () const;

        /**
         * Return a reference to the DoFHandler that is used to discretize
         * the temperature equation.
         */
        const parallel::distributed::DoFhandler<dim> &
        get_temperature_dof_handler () const;
        /** @} */
    };
    
  }
}


#endif
