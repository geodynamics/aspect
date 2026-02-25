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

#include <aspect/mesh_deformation/interface.h>
#include <aspect/mesh_deformation/external_tool_interface.h>
#include <aspect/simulator_access.h>

#include <aspect/python_helper.h>

namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * A plugin based on the ExternalToolInterface that loads
     * a Python script and passes the necessary information
     * via a specific, simple interface to it. It is written
     * to support coupling with the landscape evaluation code
     * Landlab, but that is not strictly required as long as
     * the same functions are made available in the file loaded.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class Landlab : public ExternalToolInterface<dim>
    {
      public:
        Landlab();

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
         * How many MPI ranks are participating in the Landlab simulation?
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
         * The argument to pass to the Landlab Python module.
         */
        std::string script_argument;

        /**
         * Whether the ASPECT geometry is spherical.
         */
        bool is_spherical;
#ifdef ASPECT_WITH_PYTHON
        /**
         * The Python module object.
         */
        PyObject *pModule;
#endif
    };
  }
}

#endif
