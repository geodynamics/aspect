/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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


#ifndef _aspect_time_stepping_function_h
#define _aspect_time_stepping_function_h

#include <aspect/time_stepping/interface.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace TimeStepping
  {
    using namespace dealii;

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
