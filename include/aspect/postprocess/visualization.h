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


#ifndef _aspect_postprocess_visualization_h
#define _aspect_postprocess_visualization_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/plugins.h>

#include <deal.II/base/thread_management.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/numerics/data_out.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * Compute the arithmetic average over q for each m of the variable quantities[q](m).
       */
      inline void average_quantities(std::vector<Vector<double> > &quantities)
      {
        const unsigned int N = quantities.size();
        const unsigned int M = quantities[0].size();
        for (unsigned int m=0; m<M; ++m)
          {
            double sum = 0;
            for (unsigned int q=0; q<N; ++q)
              sum += quantities[q](m);

            const double average = sum/N;
            for (unsigned int q=0; q<N; ++q)
              quantities[q](m) = average;
          }
      }

      /**
       * This class declares the public interface of visualization
       * postprocessors. Visualization postprocessors are used to compute
       * derived data, e.g. wave speeds, friction heating terms, etc, to be
       * put into graphical output files. They are plugins for the
       * aspect::Postprocess::Visualization class.
       *
       * Classes derived from this type must implement the functions that save
       * the state of the object and restore it (for checkpoint/restart
       * capabilities) as well as functions that declare and read parameters.
       * However, this class also already provides default implementations of
       * these functions that simply do nothing. This is appropriate for
       * objects that are stateless, as is commonly the case for visualization
       * postprocessors.
       *
       * Access to the data of the simulator is granted by the @p protected
       * member functions of the SimulatorAccess class, i.e., classes
       * implementing this interface will in general want to derive from both
       * this Interface class as well as from the SimulatorAccess class.
       *
       * <h3> How visualization plugins work </h3>
       *
       * There are two ways in which visualization plugins can work to get
       * data from a simulation into an output file:
       * <ul>
       * <li> Classes derived from this class can also derive from the deal.II
       * class DataPostprocessor or any of the classes like
       * DataPostprocessorScalar or DataPostprocessorVector. These classes can
       * be thought of as filters: DataOut will call a function in them for
       * every cell and this function will transform the values or gradients
       * of the solution and other information such as the location of
       * quadrature points into the desired quantity to output. A typical case
       * would be if the quantity $g(x)$ you want to output can be written as
       * a function $g(x) = G(u(x),\nabla u(x), x, ...)$ in a point-wise sense
       * where $u(x)$ is the value of the solution vector (i.e., the
       * velocities, pressure, temperature, etc) at an evaluation point. In
       * the context of this program an example would be to output the density
       * of the medium as a spatially variable function since this is a
       * quantity that for realistic media depends point-wise on the values of
       * the solution.
       *
       * Using this way of describing a visualization postprocessor will yield
       * a class  that would then have the following base classes: -
       * aspect::Postprocess::VisualizationPostprocessors::Interface -
       * aspect::SimulatorAccess - dealii::DataPostprocessor or any of the
       * other ones listed above
       *
       * <li> The second possibility is for a class to not derive from
       * dealii::DataPostprocessor but instead from the CellDataVectorCreator
       * class. In this case, a visualization postprocessor would generate and
       * return a vector that consists of one element per cell. The intent of
       * this option is to output quantities that are not point-wise functions
       * of the solution but instead can only be computed as integrals or
       * other functionals on a per-cell basis. A typical case would be error
       * estimators that do depend on the solution but not in a point-wise
       * sense; rather, they yield one value per cell of the mesh. See the
       * documentation of the CellDataVectorCreator class for more
       * information.
       * </ul>
       * @ingroup Postprocessing
       * @ingroup Visualization
       */
      template <int dim>
      class Interface
      {
        public:
          /**
           * Destructor. Does nothing but is virtual so that derived classes
           * destructors are also virtual.
           */
          virtual
          ~Interface ();

          /**
           * Initialize function.
           */
          virtual void initialize ();

          /**
           * Update any temporary information needed by the visualization postprocessor.
           */
          virtual void update();

          /**
           * Declare the parameters this class takes through input files.
           * Derived classes should overload this function if they actually do
           * take parameters; this class declares a fall-back function that
           * does nothing, so that postprocessor classes that do not take any
           * parameters do not have to do anything at all.
           *
           * This function is static (and needs to be static in derived
           * classes) so that it can be called without creating actual objects
           * (because declaring parameters happens before we read the input
           * file and thus at a time when we don't even know yet which
           * postprocessor objects we need).
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           * The default implementation in this class does nothing, so that
           * derived classes that do not need any parameters do not need to
           * implement it.
           */
          virtual
          void
          parse_parameters (ParameterHandler &prm);

          /**
           * A function that is used to indicate to the postprocessor manager which
           * other postprocessor(s) the current one depends upon. The returned
           * list contains the names (as strings, as you would write them in
           * the input file) of the postprocessors it requires. The manager
           * will ensure that these postprocessors are indeed used, even if
           * they were not explicitly listed in the input file, and are indeed
           * run <i>before</i> this postprocessor every time they are executed.
           *
           * The postprocessors you can nominate here are of the general
           * postprocessor class, not visualization postprocessors.
           *
           * The default implementation of this function returns an empty list.
           */
          virtual
          std::list<std::string>
          required_other_postprocessors () const;


          /**
           * Save the state of this object to the argument given to this
           * function. This function is in support of checkpoint/restart
           * functionality.
           *
           * Derived classes can implement this function and should store
           * their state in a string that is deposited under a key in the map
           * through which the respective class can later find the status
           * again when the program is restarted. A legitimate key to store
           * data under is <code>typeid(*this).name()</code>. It is up to
           * derived classes to decide how they want to encode their state.
           *
           * The default implementation of this function does nothing, i.e.,
           * it represents a stateless object for which nothing needs to be
           * stored at checkpoint time and nothing needs to be restored at
           * restart time.
           *
           * @param[in,out] status_strings The object into which
           * implementations in derived classes can place their status under a
           * key that they can use to retrieve the data.
           */
          virtual
          void save (std::map<std::string, std::string> &status_strings) const;

          /**
           * Restore the state of the object by looking up a description of
           * the state in the passed argument under the same key under which
           * it was previously stored.
           *
           * The default implementation does nothing.
           *
           * @param[in] status_strings The object from which the status will
           * be restored by looking up the value for a key specific to this
           * derived class.
           */
          virtual
          void load (const std::map<std::string, std::string> &status_strings);
      };



      /**
       * As explained in the documentation of the Interface class, the second
       * kind of visualization plugin is one that wants to generate cell-wise
       * data. Classes derived from this class need to implement a function
       * execute() that computes these cell-wise values and return a pair of
       * values where the first one indicates the name of a variable and the
       * second one is a vector with one entry per cell. This class is the
       * interface that such plugins have to implement.
       *
       * @ingroup Postprocessing
       * @ingroup Visualization
       */
      template <int dim>
      class CellDataVectorCreator : public Interface<dim>
      {
        public:
          /**
           * Destructor.
           */
          ~CellDataVectorCreator ()  override = default;

          /**
           * The function classes have to implement that want to output
           * cell-wise data.
           *
           * @return A pair of values with the following meaning:
           * - The first element provides the name by which this data should
           *   be written to the output file.
           * - The second element is a pointer to a vector with one element
           *   per active cell on the current processor. Elements corresponding
           *   to active cells that are either artificial or ghost cells (in
           *   deal.II language, see the deal.II glossary)
           *   will be ignored but must nevertheless exist in the returned
           *   vector.
           * While implementations of this function must create this
           * vector, ownership is taken over by the caller of this function
           * and the caller will take care of destroying the vector pointed
           * to.
           */
          virtual
          std::pair<std::string, Vector<float> *>
          execute () const = 0;
      };


      /**
       * This class is a tag class: If a visualization postprocessor is derived
       * from it, then this is interpreted as saying that the class will only
       * be used to generate graphical output on the surface of the model,
       * rather than for the entire domain.
       */
      template <int dim>
      class SurfaceOnlyVisualization
      {
        public:
          /**
           * Destructor. Made `virtual` to ensure that it is possible to
           * test whether a derived class is derived from this class via
           * a `dynamic_cast`.
           */
          virtual
          ~SurfaceOnlyVisualization () = default;
      };
    }


    /**
     * A postprocessor that generates graphical output in periodic intervals
     * or every time step. The time interval between generating graphical
     * output is obtained from the parameter file.
     *
     * While this class acts as a plugin, i.e. as a postprocessor that can be
     * registered with the postprocessing manager class
     * aspect::Postprocess::Manager, this class at the same time also acts as
     * a manager for plugins itself, namely for classes derived from the
     * VisualizationPostprocessors::Interface class that are used to output
     * different aspects of the solution, such as for example a computed
     * seismic wave speed from temperature, velocity and pressure.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class Visualization : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Visualization ();

        /**
         * Generate graphical output from the current solution.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Update any temporary information needed by the visualization postprocessor.
         */
        void
        update () override;

        /**
         * A function that is used to register visualization postprocessor
         * objects in such a way that the Manager can deal with all of them
         * without having to know them by name. This allows the files in which
         * individual postprocessors are implement to register these
         * postprocessors, rather than also having to modify the Manage class
         * by adding the new postprocessor class.
         *
         * @param name The name under which this visualization postprocessor
         * is to be called in parameter files.
         * @param description A text description of what this visualization
         * plugin does and that will be listed in the documentation of the
         * parameter file.
         * @param declare_parameters_function A pointer to a function that
         * declares the parameters for this postprocessor.
         * @param factory_function A pointer to a function that creates such a
         * postprocessor object and returns a pointer to it.
         */
        static
        void
        register_visualization_postprocessor (const std::string &name,
                                              const std::string &description,
                                              void (*declare_parameters_function) (ParameterHandler &),
                                              VisualizationPostprocessors::Interface<dim> *(*factory_function) ());

        /**
         * A function that is used to indicate to the postprocessor manager which
         * other postprocessor(s) the current one depends upon.
         *
         * For the current class, we simply loop over all of the visualization
         * postprocessors and collect what they want.
         */
        std::list<std::string>
        required_other_postprocessors () const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

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
         * Return the value of the parameter @p pointwise_stress_and_strain
         * that is controlled by the parameter "Point-wise stress and strain".
         */
        bool output_pointwise_stress_and_strain() const;

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
         * Interval between the generation of graphical output. This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) at which the last graphical output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        double last_output_time;

        /**
         * Maximum number of steps between the generation of graphical output.
         * This parameter
         * is read from the input file and consequently is not part of the
         * state that needs to be saved and restored.
         */
        unsigned int maximum_timesteps_between_outputs;

        /**
         * Timestep at which the last graphical output was produced
         * Used to check for the next necessary output time.
         */
        unsigned int last_output_timestep;

        /**
         * Consecutively counted number indicating the how-manyth time we will
         * create output the next time we get to it.
         */
        unsigned int output_file_number;

        /**
         * Graphical output format.
         */
        std::string output_format;

        /**
         * VTU file output supports grouping files from several CPUs into one
         * file using MPI I/O when writing on a parallel filesystem. 0 means
         * no grouping (and no parallel I/O). 1 will generate one big file
         * containing the whole solution.
         */
        unsigned int group_files;

        /**
         * On large clusters it can be advantageous to first write the
         * output to a temporary file on a local file system and later
         * move this file to a network file system. If this variable is
         * set to a non-empty string it will be interpreted as a temporary
         * storage location.
         */
        std::string temporary_output_location;

        /**
         * deal.II offers the possibility to linearly interpolate output
         * fields of higher order elements to a finer resolution. This
         * somewhat compensates the fact that most visualization software only
         * offers linear interpolation between grid points and therefore the
         * output file is a very coarse representation of the actual solution
         * field. Activating this option increases the spatial resolution in
         * each dimension by a factor equal to the polynomial degree used for
         * the velocity finite element (usually 2).
         */
        bool interpolate_output;

        /**
         * deal.II offers the possibility to filter duplicate vertices in HDF5
         * output files. This merges the vertices of adjacent cells and
         * therefore saves disk space, but misrepresents discontinuous
         * output properties. Activating this function reduces the disk space
         * by about a factor of $2^{dim}$ for hdf5 output.
         */
        bool filter_output;

        /**
         * If true, return quantities related to stresses and strain with
         * point-wise values. Otherwise the values will be averaged on each
         * cell.
         */
        bool pointwise_stress_and_strain;

        /**
         * deal.II offers the possibility to write vtu files with higher order
         * representations of the output data. This means each cell will correctly
         * show the higher order representation of the output data instead of the
         * linear interpolation between vertices that ParaView and Visit usually show.
         * Note that activating this option is safe and recommended, but requires that
         * (i) ``Output format'' is set to ``vtu'', (ii) ``Interpolate output'' is
         * set to true, (iii) you use a sufficiently new version of Paraview
         * or Visit to read the files (Paraview version 5.5 or newer, and Visit version
         * to be determined), and (iv) you use deal.II version 9.1.0 or newer.
         */
        bool write_higher_order_output;

        /**
         * For mesh deformation computations ASPECT uses an Arbitrary-Lagrangian-
         * Eulerian formulation to handle deforming the domain, so the mesh
         * has its own velocity field.  This may be written as an output field
         * by setting output_mesh_velocity to true.
         */
        bool output_mesh_velocity;

        /**
         * File operations can potentially take a long time, blocking the
         * progress of the rest of the model run. Setting this variable to
         * 'true' moves this process into a background thread, while the
         * rest of the model continues.
         */
        bool write_in_background_thread;

        /**
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);

        /**
         * Record that the mesh changed. This helps some output writers avoid
         * writing the same mesh multiple times.
         */
        void mesh_changed_signal ();

        /**
         * A function that writes the text in the second argument to a file
         * with the name given in the first argument. The function is run on a
         * separate thread to allow computations to continue even though
         * writing data is still continuing. The function takes over ownership
         * of these arguments and deletes them at the end of its work.
         */
        static
        void writer (const std::string filename,
                     const std::string temporary_filename,
                     const std::string *file_contents);

        /**
         * A list of postprocessor objects that have been requested in the
         * parameter file.
         */
        std::list<std::unique_ptr<VisualizationPostprocessors::Interface<dim> > > postprocessors;

        /**
         * A structure that keeps some history about past output operations.
         * These variables are grouped into a structure because we need them
         * twice: Once for the cell output case (via DataOut) and once for
         * surface output (via DataOutFaces).
         */
        struct OutputHistory
        {
          /**
           * Constructor
           */
          OutputHistory ();

          /**
           * Destructor. Makes sure that any background thread that may still be
           * running writing data to disk finishes before the current object is
           * fully destroyed.
           */
          ~OutputHistory ();

          /**
           * Serialize the contents of this class as far as they are not read
           * from input parameter files.
           */
          template <class Archive>
          void serialize (Archive &ar, const unsigned int version);

          /**
           * Whether the mesh changed since the last time we produced cell-based
           * output.
           */
          bool mesh_changed;

          /**
           * The most recent name of the mesh file, used to avoid redundant mesh
           * output.
           */
          std::string last_mesh_file_name;

          /**
          * A list of pairs (time, pvtu_filename) that have so far been written
          * and that we will pass to DataOutInterface::write_pvd_record to
          * create a master file that can make the association between
          * simulation time and corresponding file name (this is done because
          * there is no way to store the simulation time inside the .pvtu or
          * .vtu files).
          */
          std::vector<std::pair<double,std::string> > times_and_pvtu_names;

          /**
           * A list of list of filenames, sorted by timestep, that correspond to
           * what has been created as output. This is used to create a master
           * .visit file for the entire simulation.
           */
          std::vector<std::vector<std::string> > output_file_names_by_timestep;

          /**
           * A set of data related to XDMF file sections describing the HDF5
           * heavy data files created. These contain things such as the
           * dimensions and names of data written at all steps during the
           * simulation.
           */
          std::vector<XDMFEntry>  xdmf_entries;

          /**
           * Handle to a thread that is used to write data in the background.
           * The writer() function runs on this background thread when outputting
           * data for the `data_out` object.
           */
          Threads::Thread<void> background_thread;
        };

        /**
         * Information about the history of writing graphical
         * output for cells (via DataOut).
         */
        OutputHistory cell_output_history;

        /**
         * Information about the history of writing graphical
         * output for faces (via DataOutFaces).
         */
        OutputHistory face_output_history;

        /**
         * Write the various master record files. The master files are used by
         * visualization programs to identify which of the output files in a
         * directory, possibly one file written by each processor, belong to a
         * single time step and/or form the different time steps of a
         * simulation. For Paraview, this is a <code>.pvtu</code> file per
         * time step and a <code>.pvd</code> for all time steps. For Visit it
         * is a <code>.visit</code> file per time step and one for all time
         * steps.
         *
         * @param data_out The DataOut object that was used to write the
         * solutions.
         * @param solution_file_prefix The stem of the filename to be written.
         * @param filenames List of filenames for the current output from all
         * processors.
         * @param output_history The OutputHistory object to fill.
         */
        template <typename DataOutType>
        void write_master_files (const DataOutType &data_out,
                                 const std::string &solution_file_prefix,
                                 const std::vector<std::string> &filenames,
                                 OutputHistory                  &output_history) const;


        /**
         * A function, called from the execute() function, that takes a
         * DataOut object, does some preliminary work with it, and then
         * writes the result out to files via the writer() function (in the
         * case of VTU output) or through the XDMF facilities.
         *
         * The function returns the base name of the output files produced,
         * which can then be used for the statistics file and screen output.
         */
        template <typename DataOutType>
        std::string write_data_out_data(DataOutType   &data_out,
                                        OutputHistory &output_history) const;
    };
  }


  /**
   * Given a class name, a name, and a description for the parameter file for
   * a postprocessor, register it with the aspect::Postprocess::Manager class.
   *
   * @ingroup Postprocessing
   */
#define ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::Postprocess::VisualizationPostprocessors::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::Postprocess::Visualization<2>::register_visualization_postprocessor, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::Postprocess::VisualizationPostprocessors::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::Postprocess::Visualization<3>::register_visualization_postprocessor, \
                                name, description); \
  }
}


#endif
