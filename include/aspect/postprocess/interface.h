#ifndef __aspect__postprocess_interface_h
#define __aspect__postprocess_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>
#include <aspect/material_model/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/adiabatic_conditions.h>

#include <deal.II/base/std_cxx1x/shared_ptr.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/base/conditional_ostream.h>

#include <boost/serialization/split_member.hpp>


namespace aspect
{
  using namespace dealii;

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
         * Destructor. Does nothing but is virtual so that derived classes
         * destructors are also virtual.
         **/
        virtual
        ~Interface ();

        /**
         * Execute this postprocessor. Derived classes will implement this function
         * to do whatever they want to do to evaluate the solution at the current
         * time step.
         *
         * @param[in,out] statistics An object that contains statistics that are collected
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

        /**
         * Declare the parameters this class takes through input files.
         * Derived classes should overload this function if they actually
         * do take parameters; this class declares a fall-back function
         * that does nothing, so that postprocessor classes that do not
         * take any parameters do not have to do anything at all.
         *
         * This function is static (and needs to be static in derived
         * classes) so that it can be called without creating actual
         * objects (because declaring parameters happens before we read
         * the input file and thus at a time when we don't even know yet
         * which postprocessor objects we need).
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file. The default implementation in this class does nothing,
         * so that derived classes that do not need any parameters do
         * not need to implement it.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


        /**
         * Save the state of this object to the argument given to this function.
         * This function is in support of checkpoint/restart functionality.
         *
         * Derived classes can implement this function and should store their
         * state in a string that is deposited under a key in the map through
         * which the respective class can later find the status again when the
         * program is restarted. A legitimate key to store data under is
         * <code>typeid(*this).name()</code>. It is up to derived classes to
         * decide how they want to encode their state.
         *
         * The default implementation of this function does nothing, i.e., it
         * represents a stateless object for which nothing needs to be stored
         * at checkpoint time and nothing needs to be restored at restart time.
         *
         * @param[in,out] status_strings The object into which implementations
         * in derived classes can place their status under a key that they can
         * use to retrieve the data.
         **/
        virtual
        void save (std::map<std::string, std::string> &status_strings) const;

        /**
         * Restore the state of the object by looking up a description of the
         * state in the passed argument under the same key under which it
         * was previously stored.
         *
         * The default implementation does nothing.
         *
         * @param[in] status_strings The object from which the status will
         * be restored by looking up the value for a key specific to this
         * derived class.
         **/
        virtual
        void load (const std::map<std::string, std::string> &status_strings);
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
         * Initialize this class for a given simulator.
         *
         * @param simulator A reference to the main simulator object to which the
         * postprocessor implemented in the derived class should be applied.
         **/
        void initialize (const Simulator<dim> &simulator);

      protected:
        /** @name Accessing variables that identify overall properties of the simulator */
        /** @{ */

        /**
         * Return the current simulation time in seconds.
         */
        MPI_Comm
        get_mpi_communicator () const;

        /**
         * Return a reference to the stream object that only outputs something on one
        * processor in a parallel program and simply ignores output put into it on
        * all other processors.
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
         * Return a reference to the mapping used to describe the boundary
        * of the domain.
         */
        const Mapping<dim> &
        get_mapping () const;

        /**
         * Return the directory specified in the input parameter file to be
         * the place where output files are to be placed. The string is
         * terminated by a directory separator (i.e., '/').
         */
        std::string
        get_output_directory () const;

        /**
        * Return whether things like velocities should be converted from
        * the seconds in the MKS system to years. The value of this flag
        * is set by the corresponding entry in the input parameter file.
        */
        bool
        convert_output_to_years () const;

        /**
        * Compute the error indicators in the same way they are normally used
        * for mesh refinement. The mesh is not refined when doing so, but the
        * indicators can be used when generating graphical output to check
        * why mesh refinement is proceeding as it is.
        */
        void
        get_refinement_criteria(Vector<float> &estimated_error_per_cell) const;
        /** @} */


        /** @name Accessing variables that identify the solution of the problem */
        /** @{ */


        /**
         * Return a reference to the vector that has the current
         * solution of the entire system, i.e. the velocity and
         * pressure variables as well as the temperature.  This vector
         * is associated with the DoFHandler object returned by
         * get_dof_handler().
         *
         * @note In general the vector is a distributed vector; however, it
         * contains ghost elements for all locally relevant degrees of freedom.
         */
        const LinearAlgebra::BlockVector &
        get_solution () const;

        /**
         * Return a reference to the vector that has the solution
         * of the entire system at the previous time step.
         * This vector is associated with the DoFHandler object returned by
         * get_stokes_dof_handler().
         *
         * @note In general the vector is a distributed vector; however, it
         * contains ghost elements for all locally relevant degrees of freedom.
         */
        const LinearAlgebra::BlockVector &
        get_old_solution () const;

        /**
         * Return a reference to the DoFHandler that is used to discretize
         * the variables at the current time step.
         */
        const DoFHandler<dim> &
        get_dof_handler () const;

        /**
         * Fill the argument with a set of depth averages of the current
         * temperature field. The function fills a vector that contains
         * average temperatures over slices of the domain of same depth. The
         * function resizes the output vector to match the number of depth
         * slices.
         */
        void
        get_depth_average_temperature(std::vector<double> &values) const;

        /**
         * Compute a depth average of the current viscosity
         */
        void
        get_depth_average_viscosity(std::vector<double> &values) const;

        /**
         * Compute a depth average of the current velocity magnitude
         */
        void
        get_depth_average_velocity_magnitude(std::vector<double> &values) const;

        /**
         * Compute a depth average of the current sinking velocity
         */
        void
        get_depth_average_sinking_velocity(std::vector<double> &values) const;

        /**
         * Compute the seismic shear wave speed, Vs anomaly per element
         */
        void
        get_Vs_anomaly(Vector<float> &values) const;

        /**
         * Compute the seismic pressure wave speed, Vp anomaly per element
         */
        void
        get_Vp_anomaly(Vector<float> &values) const;

        /**
         * Compute a depth average of the seismic shear wave speed: Vs
         */
        void
        get_depth_average_Vs(std::vector<double> &values) const;
        /** @} */

        /**
         * Compute a depth average of the seismic pressure wave speed: Vp
         */
        void
        get_depth_average_Vp(std::vector<double> &values) const;
        /** @} */


        /** @name Accessing variables that identify aspects of the simulation */
        /** @{ */

        /**
         * Return a pointer to the material model to access functions like density().
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
         * Return a pointer to the object that describes the adiabatic conditions.
         */
        const AdiabaticConditions<dim> &
        get_adiabatic_conditions () const;

        /**
         * Return a pointer to the object that describes the temperature boundary
         * values.
         */
        const BoundaryTemperature::Interface<dim> &
        get_boundary_temperature () const;

        /** @} */

      private:
        /**
         * A pointer to the simulator object to which we want to
         * get access.
         */
        SmartPointer<const Simulator<dim> > simulator;
    };


    /**
     * A class that manages all objects that provide functionality to postprocess
     * solutions. It declares run time parameters for input files, reads their
     * values from such an input file, manages a list of all postprocessors selected
     * in the input file, and upon request through the execute() function calls
     * them in turn.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Manager
    {
      public:
        /**
         * Initialize the postprocessors handled by this object for a given simulator.
         *
         * @param simulator A reference to the main simulator object to which the
         * postprocessor implemented in the derived class should be applied.
         **/
        void initialize (const Simulator<dim> &simulator);

        /**
         * Execute all of the postprocessor objects that have been
         * requested in the input file. These objects also fill the
         * contents of the statistics object.
         *
         * The function returns a concatenation of the text returned by
         * the individual postprocessors.
         */
        std::list<std::pair<std::string,std::string> >
        execute (TableHandler &statistics);

        /**
         * Declare the parameters of all known postprocessors, as
         * well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file. This determines which postprocessor objects will be
         * created; then let these objects read their parameters as
         * well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Write the data of this object to a
         * stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void save (Archive &ar,
                   const unsigned int version) const;

        /**
         * Read the data of this object from a
         * stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void load (Archive &ar,
                   const unsigned int version);

        BOOST_SERIALIZATION_SPLIT_MEMBER()


        /**
         * A function that is used to register postprocessor objects
         * in such a way that the Manager can deal with all of them
         * without having to know them by name. This allows the files
         * in which individual postprocessors are implement to register
         * these postprocessors, rather than also having to modify the
         * Manage class by adding the new postprocessor class.
         *
         * @param name The name under which this postprocessor is to
         * be called in parameter files.
        * @param description A text description of what this model
        * does and that will be listed in the documentation of
        * the parameter file.
         * @param declare_parameters_function A pointer to a function
         * that declares the parameters for this postprocessor.
         * @param factory_function A pointer to a function that creates
         * such a postprocessor object and returns a pointer to it.
         **/
        static
        void
        register_postprocessor (const std::string &name,
                                const std::string &description,
                                void (*declare_parameters_function) (ParameterHandler &),
                                Interface<dim> *(*factory_function) ());

        /**
         * Exception.
         */
        DeclException1 (ExcPostprocessorNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered postprocessors.");
      private:
        /**
         * A list of postprocessor objects that have been requested
         * in the parameter file.
         */
        std::list<std_cxx1x::shared_ptr<Interface<dim> > > postprocessors;
    };


    /* -------------------------- inline and template functions ---------------------- */

    template <int dim>
    template <class Archive>
    void Manager<dim>::save (Archive &ar,
                             const unsigned int) const
    {
      // let all the postprocessors save their data in a map and then
      // serialize that
      std::map<std::string,std::string> saved_text;
      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        (*p)->save (saved_text);

      ar &saved_text;
    }


    template <int dim>
    template <class Archive>
    void Manager<dim>::load (Archive &ar,
                             const unsigned int)
    {
      // get the map back out of the stream; then let the postprocessors
      // that we currently have get their data from there. note that this
      // may not be the same set of postprocessors we had when we saved
      // their data
      std::map<std::string,std::string> saved_text;
      ar &saved_text;

      for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::iterator
           p = postprocessors.begin();
           p != postprocessors.end(); ++p)
        (*p)->load (saved_text);
    }


    /**
     * Given a class name, a name, and a description for the parameter file for a postprocessor, register it with
     * the aspect::Postprocess::Manager class.
     *
     * @ingroup Postprocessing
     */
#define ASPECT_REGISTER_POSTPROCESSOR(classname,name,description) \
  namespace ASPECT_REGISTER_POSTPROCESSOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<Interface<deal_II_dimension>,classname<deal_II_dimension> > \
    dummy_ ## classname (&aspect::Postprocess::Manager<deal_II_dimension>::register_postprocessor, \
                         name, description); }
  }
}


#endif
