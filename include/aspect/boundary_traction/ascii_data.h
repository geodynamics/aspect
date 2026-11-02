// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_traction_ascii_data_h
#define _aspect_boundary_traction_ascii_data_h

#include <aspect/boundary_traction/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace BoundaryTraction
  {
    /**
     * A class that implements prescribed traction boundary conditions from
     * data given in an AsciiData input file.
     *
     * This class supports prescribing either:
     * - Only pressure at the boundary (default behavior), where the traction is
     *   computed as the negative pressure times the outward normal vector.
     * - Full traction vector components on the boundary.
     *
     * The behavior is controlled by the input parameter
     * "Prescribe pressure instead of full traction" (default: true).
     * If true, only one component (pressure) is read from the input file and
     * used to compute traction as -pressure * normal_vector.
     * If false, the full traction vector components (equal to the spatial dimension)
     * are read from the input file.
     *
     * Additionally, the traction vector components may be specified in either
     * Cartesian coordinates or spherical coordinates, controlled by the
     * "Use spherical unit vectors" parameter.
     *
     * This flexibility allows compatibility with existing models prescribing
     * pressure boundary conditions, while enabling models requiring full traction
     * vectors.
     *
     * @ingroup BoundaryTractions
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * Return the boundary traction as a function of position. The
         * (outward) normal vector to the domain is also provided as
         * a second argument.
         */
        Tensor<1,dim>
        boundary_traction (const types::boundary_id boundary_indicator,
                           const Point<dim> &position,
                           const Tensor<1,dim> &normal_vector) const override;

        /**
         * A function that is called at the beginning of each time step to
         * indicate what the model time is for which the boundary values will
         * next be evaluated. For the current class, the function passes to
         * the parsed function what the current time is.
         */
        void
        update () override;

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
        std::set<types::boundary_id> boundary_ids;

        /**
         * Whether to specify traction in x, y, z components, or
         * r, phi, theta components.
         *
         * This variable is read from the parameter file through a parameter called 'Use spherical unit vectors'.
         */
        bool use_spherical_unit_vectors;

        /**
         * Whether to prescribe pressure (default: true) or full traction vector (false)
         * at the boundary. If true, only 1 component will be used for the boundary condition.
         *
         * This variable is read from the parameter file through a parameter called 'Prescribe pressure instead of full traction'.
         */
        bool prescribe_pressure_instead_of_full_traction;

    };
  }
}


#endif
