// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2018-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_time_stepping_function_h
#define _aspect_time_stepping_function_h

#include <aspect/time_stepping/interface.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace TimeStepping
  {
    /**
     * A class that implements a time stepping plugin on a function
     * description provided in the input file.
     *
     * @ingroup TimeStepping
     */
    template <int dim>
    class Function : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        Function () = default;

        /**
         * The main execute() function.
         */
        double
        execute() override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The function object.
         */
        Functions::ParsedFunction<1> function;
    };
  }
}


#endif
