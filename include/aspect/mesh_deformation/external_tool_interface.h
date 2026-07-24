/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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
     * parallel. This is because the points at which the external tool
     * wants to know the solution on any given MPI process will, in
     * general, not be located on the part of ASPECT's mesh that is
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
     * provides a high-level implementation of this function as part
     * of this class. In essence
     **/
    template <int dim>
    class ExternalToolInterface : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Main routine handling the mesh deformation
         */
        virtual
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_ids) const override;

        /**
         * Given all solution variables at each surface point, compute velocities at these points.
         *
         * This needs to be implemented by the derived class and implement the "surface evolution".
         */
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
         * Interpolate from velocities given in the evaluation points
         * to ASPECT velocities on the surface in form of a finite
         * element field.
         *
         * The @p velocities, which are typically vertical, were previously
         * computed by the external tool and belong to the current process.
         * The Output is a (global) finite element field vector that in
         * the velocity components of surface nodes corresponds to an
         * interpolation of the velocities provided.
         *
         * The output of this function can then be used as input for a
         * function that implements the
         * compute_velocity_constraints_on_boundary() function of the
         * base class.
         */
        LinearAlgebra::Vector
        interpolate_external_velocities_to_surface_support_points (const std::vector<Tensor<1,dim>> &velocities) const;

        /**
         * The list of evaluation points owned by the current process. These are the points where the
         * external tool will receive the ASPECT solution values and return the updated velocities.
         */
        std::vector<Point<dim>> evaluation_points;

        /**
         * deal.II RemotePointEvaluation object used to do point evaluation in the evaluation points.
         */
        std::unique_ptr<Utilities::MPI::RemotePointEvaluation<dim, dim>> remote_point_evaluator;

        /**
         * A struct to map between DoF indices and evaluation points
         */
        struct DofToEvalPointData
        {
          types::global_dof_index dof_index;
          unsigned int            evaluation_point_rank;
          unsigned int            evaluation_point_index;
          unsigned int            component;
          double                  squared_distance;

          template <class Archive>
          void
          serialize(Archive &ar, const unsigned int /*version*/)
          {
            ar &dof_index;
            ar &evaluation_point_rank;
            ar &evaluation_point_index;
            ar &component;
            ar &squared_distance;
          }
        };

        /**
         * A vector to store a map between Dof indices and evaluation points.
         *
         * The map will contain an entry for each DoF on the surface of the ASPECT mesh and contains
         * the index of the evaluation point that is closest to the DoF. As each support point has
         * several components (x,y,z velocity), the map contains one entry for each component.
         *
         * In a parallel computation, this map only contains entries for evaluation points owned by the
         * current process. Note that the DoF indices are not necessarily locally owned.
         *
         * This map is used in interpolate_external_velocities_to_surface_support_points() to copy
         * external velocities to each surface DoF from the closest evaluation point.
         */
        std::vector<DofToEvalPointData> map_dof_to_eval_point;
    };
  }
}

#endif
