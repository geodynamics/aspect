// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2013-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_boundary_composition_box_h
#define _aspect_boundary_composition_box_h

#include <aspect/boundary_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a box
     * geometry.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class Box : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns user-defined constant compositions at the
         * boundaries.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;

        /**
         * This function performs some basic sanity checks on the parameter
         * values previously read from the input file.
         */
        void initialize () override;

        /**
         * Declare the parameters this class takes through input files. This
         * class declares the inner and outer boundary compositions.
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
         * The values of the various composition variables on each of the
         * 2*dim boundaries of the box.
         *
         * This variable is read from the parameter file through parameters called 'Left composition',
         * 'Right composition', 'Bottom composition', 'Top composition',
         * 'Front composition', and 'Back composition'.
         */
        std::vector<double> composition_values[2*dim];
    };
  }
}


#endif
