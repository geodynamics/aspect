/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_deformation_interface_h
#define _aspect_mesh_deformation_interface_h

#include <aspect/plugins.h>
#include <aspect/simulator_access.h>

#include <aspect/global.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_handler.h>

#if DEAL_II_VERSION_GTE(9,1,0)
#  include <deal.II/lac/affine_constraints.h>
#else
#  include <deal.II/lac/constraint_matrix.h>
#endif

#include <deal.II/base/index_set.h>


namespace aspect
{
  using namespace dealii;

  template <int dim> class Simulator;

  /**
   * A namespace that contains everything that is related to the deformation
   * of the mesh vertices over time.
   */
  namespace MeshDeformation
  {
    /**
     * A base class for mesh deformation plugins. Each derived class should
     * implement a function that determines the deformation velocity for certain
     * mesh vertices and store them in a AffineConstraints<double> object. The velocities
     * for all non-constrained vertices will be computed by solving a Laplace-
     * problem with the given constraints.
     */
    template<int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface() = default;

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         *
         * The default implementation of this function does nothing.
         */
        virtual void initialize ();

        /**
         * A function that is called at the beginning of each time step and
         * that allows the implementation to update internal data structures.
         * This is useful, for example, if you have mesh deformation that
         * depends on time, or on the solution of the previous step.
         *
         * The default implementation of this function does nothing.
         */
        virtual void update();

        /**
         * A function that returns the initial deformation of points on the
         * boundary (e.g. the surface vertices). @p position is the undeformed
         * position and this function is expected to return the
         * displacement vector of this position. The default implementation
         * returns a zero displacement (= no initial deformation).
         */
        virtual
        Tensor<1,dim>
        compute_initial_deformation_on_boundary(const types::boundary_id boundary_indicator,
                                                const Point<dim> &position) const;

        /**
         * A function that creates constraints for the velocity of certain mesh
         * vertices (e.g. the surface vertices) for a specific set of boundaries.
         * The calling class will respect
         * these constraints when computing the new vertex positions.
         * The default implementation creates no constraints.
         */
        virtual
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const;

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
    };



    /**
     * The MeshDeformationHandler that handles the motion
     * of the surface, the internal nodes and computes the
     * Arbitrary-Lagrangian-Eulerian correction terms.
     */
    template<int dim>
    class MeshDeformationHandler: public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the mesh deformation handler, allowing it to read in
         * relevant parameters as well as giving it a reference to the
         * Simulator that owns it, since it needs to make fairly extensive
         * changes to the internals of the simulator.
         */
        MeshDeformationHandler(Simulator<dim> &simulator);

        /**
         * Destructor for the mesh deformation handler.
         */
        ~MeshDeformationHandler() override;

        /**
         * Initialization function of the MeshDeformationHandler.
         *
         * The default implementation of this function does nothing.
         */
        void initialize();

        /**
         * Update function of the MeshDeformationHandler. This function
         * allows the individual mesh deformation objects to update.
         */
        void update();

        /**
         * The main execution step for the mesh deformation implementation. This
         * computes the motion of the surface, moves the boundary nodes
         * accordingly, redistributes the internal nodes in order to
         * preserve mesh regularity, and calculates the Arbitrary-
         * Lagrangian-Eulerian correction terms for advected quantities.
         */
        void execute();

        /**
         * Allocates and sets up the members of the MeshDeformationHandler. This
         * is called by Simulator<dim>::setup_dofs()
         */
        void setup_dofs();

        /**
         * Declare parameters for the mesh deformation handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the mesh deformation handling.
         */
        void parse_parameters (ParameterHandler &prm);

        /**
         * A function that is used to register mesh deformation objects in such
         * a way that the Manager can deal with all of them without having to
         * know them by name. This allows the files in which individual
         * plugins are implemented to register these plugins, rather than also
         * having to modify the Manager class by adding the new initial
         * mesh deformation plugin class.
         *
         * @param name A string that identifies the mesh deformation model
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this mesh deformation model
         * wants to read from input files.
         * @param factory_function A pointer to a function that can create an
         * object of this mesh deformation model.
         */
        static
        void
        register_mesh_deformation
        (const std::string &name,
         const std::string &description,
         void (*declare_parameters_function) (ParameterHandler &),
         Interface<dim> *(*factory_function) ());

        /**
         * Return a map of boundary indicators to the names of all mesh deformation models currently
         * used in the computation, as specified in the input file.
         */
        const std::map<types::boundary_id, std::vector<std::string> > &
        get_active_mesh_deformation_names () const;

        /**
         * Return a map of boundary indicators to vectors of pointers to all mesh deformation models
         * currently used in the computation, as specified in the input file.
         */
        const std::map<types::boundary_id,std::vector<std::unique_ptr<Interface<dim> > > > &
        get_active_mesh_deformation_models () const;

        /**
         * Return a set of all the indicators of boundaries with
         * mesh deformation objects on them.
         */
        const std::set<types::boundary_id> &
        get_active_mesh_deformation_boundary_indicators () const;

        /**
         * Return the boundary id of the surface that has a free surface
         * mesh deformation object. If no free surface is used,
         * an empty set is returned.
         */
        const std::set<types::boundary_id> &
        get_free_surface_boundary_indicators () const;

        /**
         * Return the initial topography stored on
         * the Q1 finite element that describes the mesh geometry.
         * Note that a topography is set for all mesh nodes,
         * but only the values of surface boundary nodes are correct.
         * The internal nodes get the same initial topography as the
         * corresponding surface node. In other words, there is
         * no decrease of the initial topography with depth.
         * However, only the topography stored at the surface nodes
         * is taken into account in the diffusion plugin that
         * uses this function. TODO Once all initial_topography
         * is prescribed through initial_mesh_deformation, this
         * function can be removed.
         */
        const LinearAlgebra::Vector &
        get_initial_topography () const;

        /**
         * Return the mesh displacements stored on
         * the mesh deformation element.
         */
        const LinearAlgebra::Vector &
        get_mesh_displacements () const;

        /**
         * Go through the list of all mesh deformation objects that have been selected
         * in the input file (and are consequently currently active) and return
         * true if one of them has the desired type specified by the template
         * argument.
         */
        template <typename MeshDeformationType>
        bool
        has_matching_mesh_deformation_object () const;

        /**
         * Go through the list of all mesh deformation objects that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be casted to that type. If so, return a reference
         * to it. If no mesh deformation object is active that matches the given type,
         * throw an exception.
         */
        template <typename MeshDeformationType>
        const MeshDeformationType &
        get_matching_mesh_deformation_object () const;

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
        DeclException1 (ExcMeshDeformationNameNotFound,
                        std::string,
                        << "Could not find entry <"
                        << arg1
                        << "> among the names of registered mesh deformation objects.");

      private:
        /**
        * Set the boundary conditions for the solution of the elliptic
        * problem, which computes the initial displacements of the internal
        * vertices so that the mesh does not become too distorted due to
        * motion of the surface. Displacements of vertices on the
        * deforming surface are fixed according to the selected deformation
        * plugins.
        */
        AffineConstraints<double> make_initial_constraints ();

        /**
         * Deform the initial mesh by solving a Laplace equation
         * for the interior mesh vertices. The boundary deformation
         * is prescribed as given by the
         * compute_initial_deformation_on_boundary() function of
         * the individual mesh deformation plugins.
         */
        void deform_initial_mesh ();

        /**
         * Set the boundary conditions for the solution of the elliptic
         * problem, which computes the displacements of the internal
         * vertices so that the mesh does not become too distorted due to
         * motion of the surface. Velocities of vertices on the
         * deforming surface are fixed according to the selected deformation
         * plugins. Velocities of vertices on free-slip boundaries are
         * constrained to be tangential to those boundaries. Velocities of
         * vertices on no-slip boundaries are set to be zero. If a no-slip
         * boundary is marked as additional tangential, then vertex velocities
         * are constrained as tangential.
         */
        void make_constraints ();

        /**
         * Solve vector Laplacian equation for internal mesh displacements.
         */
        void compute_mesh_displacements ();

        /**
         * Set up the vector with initial displacements of the mesh
         * due to the initial topography, as supplied by the initial
         * topography plugin based on the surface coordinates of the
         * mesh nodes. We set all entries to the initial topography
         * based on its surface coordinates, i.e. the initial topography
         * is not corrected for depth from the surface as it is
         * for the initial mesh deformation. TODO this is ok for now,
         * because the surface diffusion plugin only cares about the
         * initial topography at the surface, but it would be more correct if it
         * sets the initial topography to the actual initial distortion of
         * the mesh cells. When all initial_topography plugins are converted
         * to the new initial_mesh_deformation functionality, this function
         * can be removed.
         */
        void set_initial_topography ();

        /**
         * Calculate the velocity of the mesh for ALE corrections.
         */
        void interpolate_mesh_velocity ();

        /**
         * Reference to the Simulator object to which a MeshDeformationHandler
         * instance belongs.
         */
        Simulator<dim> &sim;

        /**
        * Finite element for the mesh deformation implementation, which is
        * used for tracking mesh deformation.
        */
        const FESystem<dim> mesh_deformation_fe;

        /**
         * DoFHandler for the mesh deformation implementation.
         */
        DoFHandler<dim> mesh_deformation_dof_handler;

        /**
         * BlockVector which stores the mesh velocity in the
         * Stokes finite element space.
         * This is used for ALE corrections.
         */
        LinearAlgebra::BlockVector mesh_velocity;

        /**
         * Vector for storing the positions of the mesh vertices. This
         * is used for calculating the mapping from the reference cell to
         * the position of the cell in the deformed mesh. This must be
         * redistributed upon mesh refinement.
         */
        LinearAlgebra::Vector mesh_displacements;

        /**
         * Vector for storing the positions of the mesh vertices at the initial timestep.
         * This must be redistributed upon mesh refinement.
         * We need to store the initial topography because it is not taken
         * into account into the mesh displacements used by the MappingQ1Eulerian.
         * The current mesh displacements plus the initial topography provide
         * the actual topography at any time.
         */
        LinearAlgebra::Vector initial_topography;

        /**
         * Vector for storing the mesh velocity in the mesh deformation finite
         * element space, which is, in general, not the same finite element
         * space as the Stokes system. This is used for interpolating
         * the mesh velocity in the mesh deformation finite element space onto
         * the velocity in the Stokes finite element space, which is then
         * used for making the ALE correction in the advection equations.
         */
        LinearAlgebra::Vector fs_mesh_velocity;

        /**
         * IndexSet for the locally owned DoFs for the mesh system
         */
        IndexSet mesh_locally_owned;

        /**
         * IndexSet for the locally relevant DoFs for the mesh system
         */
        IndexSet mesh_locally_relevant;

        /**
         * Storage for the mesh velocity constraints for solving the
         * elliptic problem.
         */
        AffineConstraints<double> mesh_velocity_constraints;

        /**
         * Storage for the mesh vertex constraints for keeping the mesh conforming
         * upon redistribution.
         */
        AffineConstraints<double> mesh_vertex_constraints;

        /**
         * A map of boundary ids to mesh deformation objects that have been requested
         * in the parameter file.
         */
        std::map<types::boundary_id,std::vector<std::unique_ptr<Interface<dim> > > > mesh_deformation_objects;

        /**
         * Map from boundary id to a vector of names representing
         * mesh deformation objects.
         */
        std::map<types::boundary_id, std::vector<std::string> > mesh_deformation_object_names;

        /**
         * The set of boundary indicators for which mesh deformation
         * objects are set and that therefore can deform over time as
         * prescribed in the mesh_deformation_objects.
         */
        std::set<types::boundary_id> prescribed_mesh_deformation_boundary_indicators;

        /**
         * A set of boundary indicators that denote those boundaries that are
         * allowed to move their mesh tangential to the boundary.
         */
        std::set<types::boundary_id> tangential_mesh_deformation_boundary_indicators;

        /**
         * A set of boundary indicators, on which mesh deformation is prescribed to
         * be zero (fixed boundaries that never move). All boundaries except those
         * in prescribed_mesh_deformation_boundary_indicators and
         * tangential_mesh_deformation_boundary_indicators are in this set.
         */
        std::set<types::boundary_id> zero_mesh_deformation_boundary_indicators;

        /**
         * The boundary indicator(s) of the free surface(s). This is the
         * subset of prescribed_mesh_deformation_boundary_indicators for which
         * the 'free surface' plugin was selected.
         */
        std::set<types::boundary_id> free_surface_boundary_indicators;

        bool include_initial_topography;

        friend class Simulator<dim>;
        friend class SimulatorAccess<dim>;
    };



    template <int dim>
    template <typename MeshDeformationType>
    inline
    bool
    MeshDeformationHandler<dim>::has_matching_mesh_deformation_object () const
    {
      for (typename std::map<types::boundary_id, std::vector<std::unique_ptr<Interface<dim> > > >::iterator boundary_id
           = mesh_deformation_objects.begin();
           boundary_id != mesh_deformation_objects.end(); ++boundary_id)
        for (const auto &p : boundary_id->second)
          if (Plugins::plugin_type_matches<MeshDeformationType>(*p))
            return true;

      return false;
    }



    template <int dim>
    template <typename MeshDeformationType>
    inline
    const MeshDeformationType &
    MeshDeformationHandler<dim>::get_matching_mesh_deformation_object () const
    {
      AssertThrow(has_matching_mesh_deformation_object<MeshDeformationType> (),
                  ExcMessage("You asked MeshDeformation::MeshDeformationHandler::get_matching_mesh_deformation_object() for a "
                             "mesh deformation object of type <" + boost::core::demangle(typeid(MeshDeformationType).name()) + "> "
                             "that could not be found in the current model. Activate this "
                             "mesh deformation in the input file."));

      for (typename std::map<types::boundary_id, std::vector<std::unique_ptr<Interface<dim> > > >::iterator boundary_id
           = mesh_deformation_objects.begin();
           boundary_id != mesh_deformation_objects.end(); ++boundary_id)
        {
          typename std::vector<std::unique_ptr<Interface<dim> > >::const_iterator mesh_def;
          for (const auto &p : boundary_id->second)
            {
              if (Plugins::plugin_type_matches<MeshDeformationType>(*p))
                return Plugins::get_plugin_as_type<MeshDeformationType>(*p);
              else
                // We will never get here, because we had the Assert above. Just to avoid warnings.
                return Plugins::get_plugin_as_type<MeshDeformationType>(*(*mesh_def));
            }
        }
    }


    /**
     * Return a string that consists of the names of mesh deformation models that can
     * be selected. These names are separated by a vertical line '|' so
     * that the string can be an input to the deal.II classes
     * Patterns::Selection or Patterns::MultipleSelection.
     */
    template <int dim>
    std::string
    get_valid_model_names_pattern ();



    /**
     * Given a class name, a name, and a description for the parameter file
     * for a mesh deformation model, register it with the functions that can
     * declare their parameters and create these objects.
     *
     * @ingroup MeshDeformations
     */
#define ASPECT_REGISTER_MESH_DEFORMATION_MODEL(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_MESH_DEFORMATION_MODEL_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshDeformation::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::MeshDeformation::MeshDeformationHandler<2>::register_mesh_deformation, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshDeformation::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::MeshDeformation::MeshDeformationHandler<3>::register_mesh_deformation, \
                                name, description); \
  }
  }
}

#endif
