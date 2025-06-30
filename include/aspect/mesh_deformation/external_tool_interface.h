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


#ifndef _aspect_mesh_deformation_external_tool_interface_h
#define _aspect_mesh_deformation_external_tool_interface_h

#include <aspect/mesh_deformation/interface.h>


namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * This class provides some support for writing classes derived
     * from MeshDeformation::Interface when implementing surface
     * deformation models that are based on external tools such as
     * Fastscape, Landlab, OpenLEM, etc. The primary complication of
     * interacting with these external tools is that they generally
     * use their own discretization of surface processes, using a mesh
     * that, in most cases, will be different from the one used by
     * ASPECT. (Of course, ASPECT's use is a *volume* mesh, whereas
     * surface processes use *surface* meshes. ASPECT's own surface
     * diffusion, for example, works on the surface faces of ASPECT's
     * mesh, but in general, external tools use a mesh that is
     * *entirely unrelated* to the surface cells of ASPECT's volume
     * mesh.) In these cases, one needs to implement transfer
     * operations to take ASPECT's solution and interpolate it to the
     * external tool's mesh, run some surface evolution steps, and
     * then interpolate back to the surface nodes of ASPECT's
     * mesh. These interpolation operations are difficult to implement
     * in their own right, but are made particularly cumbersome in
     * parallel because then, the points at which the external tool
     * wants to know the solution on any given MPI process will, in
     * general, not be located on that part of ASPECT's mesh that is
     * owned by that MPI process. As a consequence, implementing the
     * interpolation requires finding the MPI process that owns the
     * cell on which that point is located, and then communication
     * back and forth.
     *
     * @sect3{Helper functions}
     *
     * This class implements the interpolation functionality for
     * derived classes to use. It works through a two-stage approach:
     * First, derived classes declare at which points they require
     * ASPECT's solution to be evaluated, via the
     * set_evaluation_points() function. This function then finds
     * which process owns the cell around this point, and sets up
     * communication structures that will make later evaluation
     * efficient. Second, this class provides the
     * evaluate_aspect_variables_at_points() function that uses these
     * communication structures to evaluate ASPECT's current solution
     * at the points previously set. The class also provides the
     * interpolate_vertical_velocities_to_surface_points() function
     * that takes a set of (vertical) velocity values at these points
     * and uses them to interpolate the information back onto ASPECT's
     * mesh.
     *
     * All three of these functions are `protected` member functions
     * of this class, ready to be called by derived classes.
     *
     *
     * @sect3{A high-level function}
     *
     * All classes derived from Interface need to implement the
     * Interface::compute_velocity_constraints_on_boundary()
     * function. Because the basic outline of this function looks
     * essentially the same for all external tools, this class also
     * provides a high-level implementation of this function is also
     * provided as part of this class. In essence, it looks as
     * follows (pseudo-code), relying on the
     * `compute_updated_velocities_at_points()` function that gets
     * ASPECT's solution at the evaluation points and returns velocity
     * vectors at all of these evaluation points:
     * @code
     *   // Interpolate ASPECT's solution at evaluation points:
     *   aspect_surface_velocities = evaluate_aspect_variables_at_points();
     *
     *   // Call derive class's method to compute updated velocities:
     *   external_surface_velocities = compute_updated_velocities_at_points(aspect_surface_velocities);
     *
     *   // Turn the result into the constraints that
     *   // compute_velocity_constraints_on_boundary is supposed
     *   // to return.
     * @endcode
     */
    template <int dim>
    class ExternalToolInterface : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        virtual
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_ids) const override;

        virtual
        std::vector<Tensor<1,dim>>
        compute_updated_velocities_at_points (const std::vector<std::vector<double>> &current_solution_at_points) const = 0;


      protected:
        // Helper functions to be used in derived classes:

        /**
         * Declare at which points the external tool driven by a class
         * derived from the current one needs to evaluate the ASPECT
         * solution. Each process will call this function with the points
         * it needs. Different processes will, in general, provide
         * distinct sets of points.
         *
         * The points provided here are given in the undeformed
         * configuration that corresponds to the initial mesh. For
         * example, if the mesh describes a box geometry, then the
         * points should like in the $x$-$y$-plane that forms the top
         * surface, rather than on the currently deformed top surface
         * that describes the current elevation map.
         *
         * @note This function sets up communication structures that
         *   encode, among other things, which process owns which
         *   points.  This information becomes outdated when the
         *   ASPECT mesh changes for mesh refinement. As a
         *   consequence, the current function sets up a process that
         *   invalidates the communication structures upon mesh
         *   refinement. Derived classes therefore have to set up
         *   their own code to call set_evaluation_points() again when
         *   the ASPECT mesh changes, or perhaps simply check whether
         *   the information needs to be updated in an overload of the
         *   update() function called at the beginning of each time
         *   step.
         */
        void
        set_evaluation_points (const std::vector<Point<dim>> &evaluation_points);

        /**
         * Return the value of the ASPECT solution at the set of points
         * previously set via set_evaluation_points().
         *
         * For each point, the return object contains a vector with as
         * many components as there are ASPECT solution components.
         */
        std::vector<std::vector<double>>
        evaluate_aspect_variables_at_points () const;

        /**
         * Given velocities (typically vertical, but not always) at
         * all of the points this process has previously set via
         * set_evaluation_points(), compute a (global) finite element
         * field vector that in the velocity components of surface
         * nodes corresponds to an interpolation of the velocities
         * provided.
         *
         * The output of this function can then be used as input for a
         * function that implements the
         * compute_velocity_constraints_on_boundary() function of the
         * base class.
         */
        LinearAlgebra::Vector
        interpolate_velocities_to_surface_points (const std::vector<Tensor<1,dim>> &vertical_velocities) const;


      private:

        std::vector<Point<dim>> evaluation_points;

        // Timo: Replace by whatever type you need here
        class X {};
        std::unique_ptr<X> remote_point_evaluator;
    };
  }
}

#endif
