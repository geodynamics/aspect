/*
  Copyright (C) 2018 by the authors of the ASPECT code.

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


#ifndef _aspect_world_builder_continental_plate_h
#define _aspect_world_builder_continental_plate_h

#include <aspect/world_builder/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace WorldBuilder
  {
    using namespace dealii;

    /**
     * A class that implements a continental plate.
     *
     * @ingroup ContinentalPlate
     */
    template <int dim>
    class ContinentalPlate : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        /**
         * Returns the temperature of an object.
         */
        double temperature (const Utilities::NaturalCoordinate<3> &position,
                            const unsigned int object_number,
                            const double depth,
                            double temperature) const;

        /**
         * Returns the composition of an object.
         */
        double composition (const Utilities::NaturalCoordinate<3> &position,
                            const unsigned int n_comp,
                            const unsigned int object_number,
                            const double depth,
                            double composition) const;

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
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        std::vector<std::vector<dealii::Point<2, double> > > object_coordinates;
        // TODO: This can be done probably more memory efficient by using hashmaps,
        // but might not be faster.
        std::vector<std::string> temperature_submodules;
        std::vector<std::string> composition_submodules;
        std::vector<double> temperatures;
        std::vector<double> compositions;
        std::vector<double> temperature_depths;
        std::vector<double> composition_depths;

    };
  }
}


#endif
