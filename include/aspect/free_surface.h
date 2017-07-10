/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_free_surface_h
#define _aspect_free_surface_h

#include <aspect/simulator.h>

namespace aspect
{
  using namespace dealii;

  template <int dim>
  class FreeSurfaceHandler
  {
    public:
      /**
       * Initialize the free surface handler, allowing it to read in
       * relevant parameters as well as giving it a reference to the
       * Simulator that owns it, since it needs to make fairly extensive
       * changes to the internals of the simulator.
       */
      FreeSurfaceHandler(Simulator<dim> &, ParameterHandler &prm);

      /**
       * Destructor for the free surface handler.
       */
      ~FreeSurfaceHandler();

      /**
       * The main execution step for the free surface implementation. This
       * computes the motion of the free surface, moves the boundary nodes
       * accordingly, redistributes the internal nodes in order to
       * preserve mesh regularity, and calculates the Arbitrary-
       * Lagrangian-Eulerian correction terms for advected quantities.
       */
      void execute();

      /**
       * Allocates and sets up the members of the FreeSurfaceHandler. This
       * is called by Simulator<dim>::setup_dofs()
       */
      void setup_dofs();

      /**
       * Apply stabilization to a cell of the system matrix.  The
       * stabilization is only added to cells on a free surface.  The
       * scheme is based on that of Kaus et. al., 2010.  Called during
       * assembly of the system matrix.
       */
      void apply_stabilization (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                internal::Assembly::Scratch::StokesSystem<dim>       &scratch,
                                internal::Assembly::CopyData::StokesSystem<dim>      &data);

      /**
       * Declare parameters for the free surface handling.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Parse parameters for the free surface handling.
       */
      void parse_parameters (ParameterHandler &prm);

    private:
      /**
       * Set the boundary conditions for the solution of the elliptic
       * problem, which computes the displacements of the internal
       * vertices so that the mesh does not become too distorted due to
       * motion of the free surface.  Velocities of vertices on the free
       * surface are set to be the normal of the Stokes velocity solution
       * projected onto that surface.  Velocities of vertices on free-slip
       * boundaries are constrained to be tangential to those boundaries.
       * Velocities of vertices on no-slip boundaries are set to be zero.
       */
      void make_constraints ();

      /**
       * Project the Stokes velocity solution onto the
       * free surface. Called by make_constraints()
       */
      void project_velocity_onto_boundary (LinearAlgebra::Vector &output);

      /**
       * Solve vector Laplacian equation for internal mesh displacements.
       */
      void compute_mesh_displacements ();

      /**
       * Calculate the velocity of the mesh for ALE corrections.
       */
      void interpolate_mesh_velocity ();

      /**
       * Reference to the Simulator object to which a FreeSurfaceHandler
       * instance belongs.
       */
      Simulator<dim> &sim;

      /**
      * Finite element for the free surface implementation, which is
      * used for tracking mesh deformation.
      */
      const FESystem<dim> free_surface_fe;

      /**
       * DoFHandler for the free surface implementation.
       */
      DoFHandler<dim> free_surface_dof_handler;


      /**
       * Stabilization parameter for the free surface.  Should be between
       * zero and one. A value of zero means no stabilization.  See Kaus
       * et. al. 2010 for more details.
       */
      double free_surface_theta;

      /**
       * BlockVector which stores the mesh velocity.
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
       * Vector for storing the mesh velocity in the free surface finite
       * element space, which is, in general, not the same finite element
       * space as the Stokes system. This is used for interpolating
       * the mesh velocity in the free surface finite element space onto
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
       * Storage for the mesh displacement constraints for solving the
       * elliptic problem
       */
      ConstraintMatrix mesh_displacement_constraints;

      /**
       * Storage for the mesh vertex constraints for keeping the mesh conforming
       * upon redistribution.
       */
      ConstraintMatrix mesh_vertex_constraints;

      /**
       * A struct for holding information about how to advect the free surface.
       */
      struct SurfaceAdvection
      {
        enum Direction { normal, vertical };
      };

      /**
       * Stores whether to advect the free surface in the normal direction
       * or the direction of the local vertical.
       */
      typename SurfaceAdvection::Direction advection_direction;


      /**
       * A set of boundary indicators that denote those boundaries that are
       * allowed to move their mesh tangential to the boundary. All
       * boundaries that have tangential material velocity boundary
       * conditions are in this set by default, but it can be extended by
       * open boundaries, boundaries with traction boundary conditions, or
       * boundaries with prescribed material velocities if requested in
       * the parameter file.
       */
      std::set<types::boundary_id> tangential_mesh_boundary_indicators;

      /**
       * A handle on the connection that connects the Stokes assembler
       * signal of the main simulator object to the apply_stabilization()
       * function. We keep track of this connection because we need to
       * break it once the current free surface handler object goes out
       * of scope.
       *
       * With the current variable, the connection is broken once the
       * scoped_connection goes out of scope, i.e., when the surrounding
       * class is destroyed.
       */
      boost::signals2::scoped_connection assembler_connection;

      friend class Simulator<dim>;
      friend class SimulatorAccess<dim>;
  };
}


#endif
