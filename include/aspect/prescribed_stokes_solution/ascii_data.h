// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2011-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_prescribed_stokes_solution_ascii_data_h
#define _aspect_prescribed_stokes_solution_ascii_data_h

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace PrescribedStokesSolution
  {
    /**
     * A class that implements a prescribed velocity field determined from
     * a AsciiData input file.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * For the current class, this function returns value from the text
         * files.
         */
        void
        stokes_solution (const Point<dim> &position, Vector<double> &value) const override;

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
    };
  }
}


#endif
