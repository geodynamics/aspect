/*
  Copyright (C) 2014 - 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/



#ifndef _aspect_mesh_refinement_minimum_refinement_function_h
#define _aspect_mesh_refinement_minimum_refinement_function_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a minimum refinement level based on a
     * functional description provided in the input file.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class MinimumRefinementFunction : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /*
         * At the beginning of each time step, update the time for the
         * ParsedFunction.
         */
        virtual
        void
        update ();

        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         *
         */
        virtual
        void
        tag_additional_cells () const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * The coordinate representation to evaluate the function. Possible
         * choices are depth, cartesian and spherical.
         */
        ::aspect::Utilities::Coordinates::CoordinateSystem coordinate_system;

        /**
         * A function object representing the minimum refinement level. The
         * function always depends on 3 variables, although in the case of the
         * 'depth' coordinate system only the first is used to evaluate the
         * function.
         */
        Functions::ParsedFunction<dim> min_refinement_level;

    };
  }
}

#endif
