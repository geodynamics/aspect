/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_initial_composition_world_builder_h
#define _aspect_initial_composition_world_builder_h

#include <aspect/global.h>

#ifdef ASPECT_WITH_WORLD_BUILDER

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>

namespace WorldBuilder
{
  class World;
}


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file through the
     * World builder.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class WorldBuilder : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

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
        /**
         * A vector that specifies for each compositional field index if the world builder
         * should be evaluated for this compositional field.
         */
        std::vector<bool> relevant_compositions;

        /**
         * A pointer to the WorldBuilder object. Keeping this pointer ensures
         * that the object doesn't go away while we still need it.
         */
        std::shared_ptr<const ::WorldBuilder::World> world_builder;
    };
  }
}

#endif
#endif
