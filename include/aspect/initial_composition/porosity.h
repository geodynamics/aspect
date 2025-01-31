/*
  Copyright (C) 2017 - 2023 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_porosity_h
#define _aspect_initial_composition_porosity_h

#include <aspect/initial_composition/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements initial conditions for the porosity field
     * by computing the equilibrium melt fraction for the given initial
     * condition and reference pressure profile. Note that this plugin only
     * works if there is a compositional field called 'porosity', and the
     * used material model implements the 'MeltFractionModel' interface.
     * All compositional fields except porosity are not changed by this plugin.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Porosity : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the plugin.
         */
        void initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const override;

      private:
        /**
         * A shared pointer to the initial temperature object
         * that ensures that the current object can continue
         * to access the initial temperature object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialTemperature::Manager<dim>> initial_temperature_manager;

        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition_manager;
    };
  }
}


#endif
