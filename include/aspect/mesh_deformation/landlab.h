/*
  Copyright (C) 2025 - 2025 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_deformation_landlab_h
#define _aspect_mesh_deformation_landlab_h

#if defined(ASPECT_WITH_PYTHON) && defined(ASPECT_WITH_LANDLAB)
#include <aspect/mesh_deformation/interface.h>
#include <aspect/mesh_deformation/parallel_unstructured_interface.h>
#include <aspect/simulator_access.h>

#include <aspect/python_helper.h>

namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * A plugin based on the ParallelUnstructuredInterface that
     * loads a Python script that runs Landlab. This class passes
     * information to it via a specific set of functions that are
     * explicitly called in the Python script.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class Landlab : public ParallelUnstructuredInterface<dim>
    {
      public:
        /**
         * Initialize function, this creates the Python interpreter and
         * and loads the Landlab Python module.
         */
        void initialize() override;

        /**
         * Update function. This sets the evaluation points and creates the
         * communication data structures.
         */
        void update() override;

        /**
         * Call into the Landlab Python module to compute the updated velocities
         * at the evaluation points.
         */
        virtual
        std::vector<Tensor<1,dim>>
        compute_updated_velocities_at_points (const std::vector<std::vector<double>> &current_solution_at_points) const override;

        /**
         * Compute the initial deformation by querying the Landlab Python module for the
         * initial topography at the evaluation points, interpolating to support points, and
         * creating corresponding constraints.
         */
        void
        compute_initial_deformation_as_constraints(const Mapping<dim> &mapping,
                                                   const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                   const types::boundary_id boundary_indicator,
                                                   AffineConstraints<double> &constraints) const override;

        /**
         * Declare parameters.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The number of MPI ranks are participating in the Landlab simulation.
         * Currently this has to be 1.
         */
        unsigned int n_landlab_ranks;

        /**
         * The MPI communicator for the Landlab simulation.
         */
        MPI_Comm landlab_communicator;

        /**
         * Whether this rank is one of the ranks that runs the Landlab simulation.
         */
        bool this_rank_runs_landlab;

        /**
         * The path to the Landlab Python module.
         */
        std::string script_path;

        /**
         * The name of the Landlab Python module without the .py extension.
         */
        std::string script_module_name;

        /**
         * Whether the ASPECT geometry is spherical.
         */
        bool is_spherical;

        /**
         * The Python module object.
         */
        PyObject *pModule;
    };
  }
}

#endif
#endif
