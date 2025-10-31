/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.templates.h>
#include <aspect/simulator/assemblers/interface.h>

namespace aspect
{
  namespace Assemblers
  {
    /**
     * Apply stabilization to a cell of the system matrix. The
     * stabilization is only added to cells on a free surface. The
     * scheme is based on that of Kaus et. al., 2010. Called during
     * assembly of the system matrix.
     */
    template <int dim>
    class ApplyStabilization: public Assemblers::Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        ApplyStabilization(const double stabilization_theta);

        void
        execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch,
                 internal::Assembly::CopyData::CopyDataBase<dim> &data) const override;

      private:
        /**
         * Stabilization parameter for the free surface. Should be between
         * zero and one. A value of zero means no stabilization. See Kaus
         * et. al. 2010 for more details.
         */
        const double free_surface_theta;
    };
  }

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
     * for all non-constrained vertices will be computed by solving a Laplace
     * problem with the given constraints.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class Interface : public Plugins::InterfaceBase
    {
      public:
        /**
         * A function that will be called to check whether stabilization is needed.
         */
        virtual bool needs_surface_stabilization() const;


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
                                                 const std::set<types::boundary_id> &boundary_ids) const;
    };



    /**
     * The MeshDeformationHandler that handles the motion
     * of the surface, the internal nodes and computes the
     * Arbitrary-Lagrangian-Eulerian correction terms.
     */
    template <int dim>
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
         * Called by Simulator::set_assemblers() to allow the FreeSurface plugin
         * to register its assembler.
         */
        void set_assemblers(const SimulatorAccess<dim> &simulator_access,
                            aspect::Assemblers::Manager<dim> &assemblers) const;

        /**
         * Update function of the MeshDeformationHandler. This function
         * allows the individual mesh deformation objects to update.
         */
        void update();

        /**
         * The main execution step for the mesh deformation implementation. This
         * computes the motion of the surface, moves the boundary nodes
         * accordingly, redistributes the internal nodes in order to
         * preserve mesh regularity, and calculates the
         * Arbitrary-Lagrangian-Eulerian correction terms for advected quantities.
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
         * Write the data of this object to a stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void save (Archive &ar,
                   const unsigned int version) const;

        /**
         * Read the data of this object from a stream for the purpose of
         * serialization.
         */
        template <class Archive>
        void load (Archive &ar,
                   const unsigned int version);

        BOOST_SERIALIZATION_SPLIT_MEMBER()

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
         std::unique_ptr<Interface<dim>> (*factory_function) ());

        /**
         * Return a map of boundary indicators to the names of all mesh deformation models currently
         * used in the computation, as specified in the input file.
         */
        const std::map<types::boundary_id, std::vector<std::string>> &
        get_active_mesh_deformation_names () const;

        /**
         * Return a map of boundary indicators to vectors of pointers to all mesh deformation models
         * currently used in the computation, as specified in the input file.
         */
        const std::map<types::boundary_id,std::vector<std::unique_ptr<Interface<dim>>>> &
        get_active_mesh_deformation_models () const;

        /**
         * Return a set of all the indicators of boundaries with
         * mesh deformation objects on them.
         */
        const std::set<types::boundary_id> &
        get_active_mesh_deformation_boundary_indicators () const;

        /**
         * Return a set of all the indicators of boundaries that
         * require surface stabilization.
         */
        const std::set<types::boundary_id> &
        get_boundary_indicators_requiring_stabilization () const;

        /**
         * Return the boundary id of the surface that has a free surface
         * mesh deformation object. If no free surface is used,
         * an empty set is returned.
         */
        const std::set<types::boundary_id> &
        get_free_surface_boundary_indicators () const;

        /**
         * Return the stabilization parameter for the free surface.
         */
        double get_free_surface_theta () const;

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
         * Return the mesh displacements based on the mesh deformation
         * DoFHandler you can access with get_mesh_deformation_dof_handler().
         */
        const LinearAlgebra::Vector &
        get_mesh_displacements () const;

        /**
         * Return the DoFHandler used to represent the mesh deformation space.
         */
        const DoFHandler<dim> &
        get_mesh_deformation_dof_handler () const;

        /**
         * Go through the list of all mesh deformation objects that have been selected
         * in the input file (and are consequently currently active) and return
         * true if one of them has the desired type specified by the template
         * argument.
         *
         * This function can only be called if the given template type (the first template
         * argument) is a class derived from the Interface class in this namespace.
         */
        template <typename MeshDeformationType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,MeshDeformationType>::value>>
        bool
        has_matching_mesh_deformation_object () const;

        /**
         * Go through the list of all mesh deformation objects that have been selected
         * in the input file (and are consequently currently active) and see
         * if one of them has the type specified by the template
         * argument or can be cast to that type. If so, return a reference
         * to it. If no mesh deformation object is active that matches the given type,
         * throw an exception.
         *
         * This function can only be called if the given template type (the first template
         * argument) is a class derived from the Interface class in this namespace.
         */
        template <typename MeshDeformationType,
                  typename = typename std::enable_if_t<std::is_base_of<Interface<dim>,MeshDeformationType>::value>>
        const MeshDeformationType &
        get_matching_mesh_deformation_object () const;


        /**
         * If multilevel solvers are used, we need a mapping on each multigrid level. These
         * are automatically updated by this handler class and can be accessed with this
         * method.
         */
        const Mapping<dim> &
        get_level_mapping(const unsigned int level) const;

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
         * Compute the initial constraints for the mesh displacement
         * on the boundaries of the domain.  This is used on the mesh
         * deformation boundaries to describe a displacement (initial
         * topography) to be used during the simulation. The
         * displacement is given by the active deformation plugins.
         */
        void make_initial_constraints ();

        /**
         * Compute the constraints for the mesh velocity on the
         * boundaries of the domain.  On the mesh deformation
         * boundaries, the velocity is given by the active deformation
         * plugins.
         *
         * Velocities on free-slip boundaries are constrained to be
         * tangential to those boundaries. Velocities on no-slip
         * boundaries are set to be zero. If a no-slip boundary is
         * marked as additional tangential, then velocities are
         * constrained as tangential.
         */
        void make_constraints ();

        /**
         * Solve vector Laplacian equation for internal mesh displacements and update
         * the current displacement vector based on the solution.
         */
        void compute_mesh_displacements ();

        /**
         * Solve vector Laplacian equation using GMG for internal mesh displacements and update
         * the current displacement vector based on the solution.
         */
        void compute_mesh_displacements_gmg ();

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
         * Update the mesh deformation for the multigrid levels.
         */
        void update_multilevel_deformation ();

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
         * mesh_displacements from the last time step.
         */
        LinearAlgebra::Vector old_mesh_displacements;

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
        std::map<types::boundary_id,std::vector<std::unique_ptr<Interface<dim>>>> mesh_deformation_objects;

        /**
         * Map from boundary id to a vector of names representing
         * mesh deformation objects.
         */
        std::map<types::boundary_id, std::vector<std::string>> mesh_deformation_object_names;

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

        /**
         * The set of boundary indicators for which the mesh deformation
         * objects need surface stabilization.
         */
        std::set<types::boundary_id> boundary_indicators_requiring_stabilization;

        bool include_initial_topography;

        /**
         * Stabilization parameter for the free surface. Should be between
         * zero and one. A value of zero means no stabilization.  See Kaus
         * et. al. 2010 for more details.
         */
        double surface_theta;

        /**
         * If required, store a mapping for each multigrid level.
         */
        MGLevelObject<std::unique_ptr<Mapping<dim>>> level_mappings;

        /**
         * One vector on each multigrid level for the mesh displacement used in the mapping.
         */
        MGLevelObject<dealii::LinearAlgebra::distributed::Vector<double>> level_displacements;

        /**
         * Multigrid transfer operator for the displacements
         */
        MGTransferMF<dim, double> mg_transfer;

        /**
         * Multigrid level constraints for the displacements
         */
        MGConstrainedDoFs mg_constrained_dofs;

        friend class Simulator<dim>;
        friend class SimulatorAccess<dim>;
    };


    template <int dim>
    template <class Archive>
    void MeshDeformationHandler<dim>::save (Archive &ar,
                                            const unsigned int) const
    {
      // let all the mesh deformation plugins save their data in a map and then
      // serialize that
      //TODO: for now we assume the same plugins are active before and after restart.
      std::map<std::string,std::string> saved_text;
      for (const auto &boundary_id_and_mesh_deformation_objects : mesh_deformation_objects)
        for (const auto &p : boundary_id_and_mesh_deformation_objects.second)
          p->save (saved_text);

      ar &saved_text;
    }


    template <int dim>
    template <class Archive>
    void MeshDeformationHandler<dim>::load (Archive &ar,
                                            const unsigned int)
    {
      // get the map back out of the stream; then let the mesh deformation plugins
      // that we currently have get their data from there. note that this
      // may not be the same set ofmesh deformation plugins we had when we saved
      // their data
      std::map<std::string,std::string> saved_text;
      ar &saved_text;

      for (const auto &boundary_id_and_mesh_deformation_objects : mesh_deformation_objects)
        for (const auto &p : boundary_id_and_mesh_deformation_objects.second)
          p->load (saved_text);
    }




    template <int dim>
    template <typename MeshDeformationType, typename>
    inline
    bool
    MeshDeformationHandler<dim>::has_matching_mesh_deformation_object () const
    {
      for (const auto &object_iterator : mesh_deformation_objects)
        for (const auto &p : object_iterator.second)
          if (Plugins::plugin_type_matches<MeshDeformationType>(*p))
            return true;

      return false;
    }



    template <int dim>
    template <typename MeshDeformationType, typename>
    inline
    const MeshDeformationType &
    MeshDeformationHandler<dim>::get_matching_mesh_deformation_object () const
    {
      AssertThrow(has_matching_mesh_deformation_object<MeshDeformationType> (),
                  ExcMessage("You asked MeshDeformation::MeshDeformationHandler::get_matching_mesh_deformation_object() for a "
                             "mesh deformation object of type <" + boost::core::demangle(typeid(MeshDeformationType).name()) + "> "
                             "that could not be found in the current model. Activate this "
                             "mesh deformation in the input file."));

      for (const auto &object_iterator : mesh_deformation_objects)
        for (const auto &p : object_iterator.second)
          if (Plugins::plugin_type_matches<MeshDeformationType>(*p))
            return Plugins::get_plugin_as_type<MeshDeformationType>(*p);

      // We will never get here, because we had the Assert above. Just to avoid warnings.
      typename std::vector<std::unique_ptr<Interface<dim>>>::const_iterator mesh_def;
      return Plugins::get_plugin_as_type<MeshDeformationType>(*(*mesh_def));

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
    aspect::internal::Plugins::RegisterHelper<aspect::MeshDeformation::Interface<2>,classname<2>> \
    dummy_ ## classname ## _2d (&aspect::MeshDeformation::MeshDeformationHandler<2>::register_mesh_deformation, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::MeshDeformation::Interface<3>,classname<3>> \
    dummy_ ## classname ## _3d (&aspect::MeshDeformation::MeshDeformationHandler<3>::register_mesh_deformation, \
                                name, description); \
  }
  }
}

#endif
