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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/



#ifndef _aspect_world_feature_features_continental_plate_h
#define _aspect_world_feature_features_continental_plate_h

#include <aspect/world_builder/features/interface.h>
#include <aspect/world_builder/world.h>


namespace aspect
{
  namespace WorldBuilder
  {
    namespace Features
    {

      class ContinentalPlate : public Interface
      {
        public:
          /**
           * constructor
           */
          ContinentalPlate(WorldBuilder::World* world);

          /**
           * Destructor
           */
          ~ContinentalPlate();

          /**
           * Read in the world builder file
           */
          virtual
          void read(ptree &property_tree);

          /**
           * Returns a temperature based on the given position
           */
          virtual
          double temperature(const std::array<double,3> position,
                             const double depth,
                             const double gravity,
                             double temperature) const;



        private:
          // local parameters
          double temperature_submodule_depth;
          double temperature_submodule_temperature;
          double composition_submodule_depth;
          unsigned int composition_submodule_composition;

      };
    }
  }
}

#endif
