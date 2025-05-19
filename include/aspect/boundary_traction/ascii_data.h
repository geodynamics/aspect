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
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        std::set<types::boundary_id> boundary_ids;

        /**
         * Whether to specify traction in x, y, z components, or
         * r, phi, theta components.
         */
        bool use_spherical_unit_vectors;

        /**
        * Whether to prescribe pressure (default: true) or full traction vector (false)
        * at the boundary. If true, only 1 component will be used for the boundary condition.
        */
        bool prescribe_pressure_instead_of_full_traction;

    };
  }
}


#endif
