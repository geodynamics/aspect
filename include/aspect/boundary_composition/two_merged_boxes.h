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


#ifndef _aspect_boundary_composition_two_merged_boxes_h
#define _aspect_boundary_composition_two_merged_boxes_h

#include <aspect/boundary_composition/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace BoundaryComposition
  {
    /**
     * A class that implements a composition boundary condition for a
     * box geometry with additional boundary indicators for the
     * lithospheric part of the left and right (and front and back in 3D) boundaries.
     *
     * @ingroup BoundaryCompositions
     */
    template <int dim>
    class TwoMergedBoxes : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * This function returns constant compositions at the boundaries.
         *
         * @copydoc aspect::BoundaryComposition::Interface::boundary_composition()
         */
        double boundary_composition (const types::boundary_id boundary_indicator,
                                     const Point<dim> &position,
                                     const unsigned int compositional_field) const override;

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

        /**
         * This function performs some basic sanity checks on the parameter
         * values previously read from the input file.
         */
        void initialize () override;

      private:
        /**
         * The values of the various composition variables on each of the
         * 2*dim+2*(dim-1) boundaries of the box.
         *
         * This variable is read from the parameter file through parameters called 'Left composition',
         * 'Right composition', 'Bottom composition', 'Top composition',
         * 'Front composition', and 'Back composition', as well as 'Left composition lithosphere',
         * 'Right composition lithosphere', 'Front composition lithosphere', and 'Back composition lithosphere'.
         */
        std::vector<double> composition_values[2*dim+2*(dim-1)];
    };
  }
}


#endif
