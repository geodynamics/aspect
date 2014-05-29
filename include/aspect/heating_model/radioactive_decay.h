/*
  Copyright (C) 2014 by the authors of the ASPECT code.

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



#ifndef __aspect__heating_model_radioactive_decay_h
#define __aspect__heating_model_radioactive_decay_h

#include <aspect/heating_model/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace HeatingModel
  {
    using namespace dealii;

    /**
     * A class that implements a heating model based on radioactive decay.
     *
     * @ingroup HeatingModels
     */
    template <int dim>
    class RadioactiveDecay : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        RadioactiveDecay ();

        /**
         * Return the specific heating rate as calculated by radioactive
         * decay.
         */
        virtual
        double
        specific_heating_rate (const double,
                               const double,
                               const std::vector<double> &composition,
                               const Point<dim> &position) const;

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
         * Number of radio active heating elements.
         */
        unsigned int                   n_radio_heating_elements;

        /**
         * Store the half life of different elements.
         */
        std::vector<double>            half_decay_times;

        /**
         * Store the unit heating rate of different elements.
         */
        std::vector<double>            radioactive_heating_rates;

        /**
         * Store the initial consentration in the crust.
         */
        std::vector<double>            radioactive_initial_concentrations_crust;

        /**
         * Store the initial consentration in the mantle.
         */
        std::vector<double>            radioactive_initial_concentrations_mantle;

        /**
         * Whether crust defined by composition or depth
         */
        bool                           is_crust_defined_by_composition;

        /**
         * Depth of the crust.
         */
        double                         crust_depth;

        /**
         * Composition number of crust.
         */
        unsigned int                   crust_composition_num;
    };
  }
}


#endif
