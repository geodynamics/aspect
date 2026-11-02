// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2014-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------



#ifndef _aspect_mesh_refinement_maximum_refinement_function_h
#define _aspect_mesh_refinement_maximum_refinement_function_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a maximum refinement level based on a
     * functional description provided in the input file.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class MaximumRefinementFunction : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * At the beginning of each time step, update the time for the
         * ParsedFunction.
         */
        void
        update () override;

        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         *
         */
        void
        tag_additional_cells () const override;

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

      private:
        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         *
         * This variable is read from the parameter file through a parameter called 'Coordinate system'.
         */
        Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * A function object representing the maximum refinement level. The
         * function always depends on 3 variables, although in the case of the
         * 'depth' coordinate system only the first is used to evaluate the
         * function.
         */
        Functions::ParsedFunction<dim> max_refinement_level;

    };
  }
}

#endif
