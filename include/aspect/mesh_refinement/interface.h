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


#ifndef __aspect__mesh_refinement_interface_h
#define __aspect__mesh_refinement_interface_h

#include <aspect/global.h>
#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/std_cxx1x/shared_ptr.h>
#include <deal.II/base/table_handler.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/distributed/tria.h>


namespace aspect
{
  using namespace dealii;

  template <int dim> class Simulator;
  template <int dim> class SimulatorAccess;


  /**
   * A namespace for everything to do with the decision on how to refine the
   * mesh every few time steps.
   *
   * @ingroup MeshRefinement
   */
  namespace MeshRefinement
  {

    /**
     * This class declares the public interface of mesh refinement plugins.
     * Plugins have two different ways to influence adaptive refinement (and
     * can make use of either or both):
     *
     * First, execute() allows the plugin to specify weights for individual
     * cells that are then used to coarsen and refine (where larger numbers
     * indicate a larger error).
     *
     * Second, after cells get flagged for coarsening and refinement (using
     * the first approach), tag_additional_cells() is executed for each
     * plugin. Here the plugin is free to set or clear coarsen and refine
     * flags on any cell.
     *
     * Access to the data of the simulator is granted by the @p protected
     * member functions of the SimulatorAccess class, i.e., classes
     * implementing this interface will in general want to derive from both
     * this Interface class as well as from the SimulatorAccess class.
     *
     * @ingroup MeshRefinement
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
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual void initialize ();

        /**
         * Execute this mesh refinement criterion. The default implementation
         * sets all the error indicators to zero.
         *
         * @param[out] error_indicators A vector that for every active cell of
         * the current mesh (which may be a partition of a distributed mesh)
         * provides an error indicator. This vector will already have the
         * correct size when the function is called.
         */
        virtual
        void
        execute (Vector<float> &error_indicators) const;

        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate. The default
         * implementation does nothing.
         *
         * This function is also called during the initial global refinement
         * cycle. At this point you do not have access to solutions,
         * DoFHandlers, or finite element spaces. You can check if this is the
         * case by querying this->get_dof_handler().n_dofs() == 0.
         */
        virtual
        void
        tag_additional_cells () const;

        /**
         * Declare the parameters this class takes through input files.
         * Derived classes should overload this function if they actually do
         * take parameters; this class declares a fall-back function that does
         * nothing, so that postprocessor classes that do not take any
         * parameters do not have to do anything at all.
         *
         * This function is static (and needs to be static in derived classes)
         * so that it can be called without creating actual objects (because
         * declaring parameters happens before we read the input file and thus
         * at a time when we don't even know yet which postprocessor objects
         * we need).
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
    };






    /**
     * A class that manages all objects that provide functionality to refine
     * meshes.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Manager : public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor. Made virtual since this class has virtual member
         * functions.
         */
        virtual ~Manager ();

        /**
         * Execute all of the mesh refinement objects that have been requested
         * in the input file. The error indicators are then each individually
         * normalized and merged according to the operation specified in the
         * input file (e.g., via a plus, a maximum operation, etc).
         */
        virtual
        void
        execute (Vector<float> &error_indicators) const;

        /**
         * Apply additional refinement criteria independent of the error
         * estimate for all of the mesh refinement objects that have been
         * requested in the input file.
         */
        virtual
        void
        tag_additional_cells () const;

        /**
         * Declare the parameters of all known mesh refinement plugins, as
         * well as of ones this class has itself.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * This determines which mesh refinement objects will be created; then
         * let these objects read their parameters as well.
         */
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * A function that is used to register mesh refinement objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implement to register these plugins, rather than also
         * having to modify the Manager class by adding the new mesh
         * refinement class.
         *
         * @param name The name under which this plugin is to be called in
         * parameter files.
         * @param description A text description of what this model does and
         * that will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that
         * declares the parameters for this plugin.
         * @param factory_function A pointer to a function that creates such a
         * mesh refinement object and returns a pointer to it.
         */
        static
        void
        register_mesh_refinement_criterion (const std::string &name,
                                            const std::string &description,
                                            void (*declare_parameters_function) (ParameterHandler &),
                                            Interface<dim> *(*factory_function) ());

        /**
         * Exception.
         */
        DeclException1 (ExcMeshRefinementNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered mesh refinement objects.");
      private:
        /**
         * An enum that describes the different ways in which we can merge the
         * results of multiple mesh refinement criteria.
         */
        enum MergeOperation
        { plus, max };

        /**
         * How to merge the results of multiple mesh refinement criteria.
         */
        MergeOperation merge_operation;

        /**
         * Whether to normalize the individual refinement indicators to the
         * range $[0,1]$ before merging.
         */
        bool normalize_criteria;

        /**
         * The scaling factors that should be applied to the individual
         * refinement indicators before merging.
         */
        std::vector<double> scaling_factors;

        /**
         * A list of mesh refinement objects that have been requested in the
         * parameter file.
         */
        std::list<std_cxx1x::shared_ptr<Interface<dim> > > mesh_refinement_objects;
    };



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a mesh refinement object, register it with the
     * aspect::MeshRefinement::Manager class.
     *
     * @ingroup MeshRefinement
     */
#define ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_MESH_REFINEMENT_CRITERION_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshRefinement::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::MeshRefinement::Manager<2>::register_mesh_refinement_criterion, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshRefinement::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::MeshRefinement::Manager<3>::register_mesh_refinement_criterion, \
                                name, description); \
  }
  }
}


#endif
