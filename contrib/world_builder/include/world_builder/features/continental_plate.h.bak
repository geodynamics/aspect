/*
  Copyright (C) 2018-2024 by the authors of the World Builder code.

  This file is part of the World Builder.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef WORLD_BUILDER_FEATURES_CONTINENTAL_PLATE_H
#define WORLD_BUILDER_FEATURES_CONTINENTAL_PLATE_H


#include "world_builder/features/interface.h"
#include "world_builder/objects/surface.h"


namespace WorldBuilder
{
  class Parameters;
  class World;

  namespace Features
  {
    namespace ContinentalPlateModels
    {
      namespace Composition
      {
        class Interface;
      }  // namespace Composition
      namespace Grains
      {
        class Interface;
      }  // namespace Grains
      namespace Temperature
      {
        class Interface;
      }  // namespace Temperature
    }  // namespace ContinentalPlateModels

    /**
     * This class represents a continental plate and can implement submodules
     * for temperature and composition. These submodules determine what
     * the returned temperature or composition of the temperature and composition
     * functions of this class will be.
     */
    class ContinentalPlate final: public Interface
    {
      public:
        /**
         * constructor
         */
        ContinentalPlate(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~ContinentalPlate() override final;

        /**
         * declare and read in the world builder file into the parameters class
         */
        static
        void declare_entries(Parameters &prm,
                             const std::string &parent_name = "",
                             const std::vector<std::string> &required_entries = {});

        /**
         * Produce a JSON snippet for the schema
         */
        static
        void make_snippet(Parameters &prm);

        /**
         * declare and read in the world builder file into the parameters class
         */
        void parse_entries(Parameters &prm) override final;


        /**
         * Returns different values at a single point in one go stored in a vector of doubles.
         *
         * The properties input decides what each entry means, and the output is generated in the
         * same order as the properties input. The properties input consists of
         * a 3D array, where the first entry identifies the property and the last two entries
         * provide extra information about that property.
         *
         * Temperature is identified by 1 and no extra information is needed. So temperature
         * input usually looks like {1,0,0}. A temperature query prodoces one entry in the output
         * vector.
         *
         * Composition is identified by 2. This produces one
         * value in the output. The second entry  identifies the composition number and the third
         * number is not used. So a commposition query asking about composition 1 looks like this:
         * {2,1,0}. A composition query prodoces one entry in the output vector.
         *
         * Grains are identified by 2. The second entry is the grain composition number and the third
         * entry is the number of grains. A query about the grains, where it asks about composition 1
         * (for example enstatite) and 500 grains, looks like this: {2,1,500}.
         * A composition query prodoces n_grains*10 entries in the output vector. The first n_grains
         * entries are the sizes of all the grains, and the other 9 entries are sets of rotation
         * matrices. The rotation matrix entries are ordered [0][0],[0][1],[0][2],[1][0],[1][1],etc.
         *
         * The entries in output variable relates the index of the property to the index in the output.
         */
        void
        properties(const Point<3> &position_in_cartesian_coordinates,
                   const Objects::NaturalCoordinate &position_in_natural_coordinates,
                   const double depth,
                   const std::vector<std::array<unsigned int,3>> &properties,
                   const double gravity,
                   const std::vector<size_t> &entry_in_output,
                   std::vector<double> &output) const override final;

      private:
        /**
         * A vector containing all the pointers to the temperature models. This vector is
         * responsible for the features and has ownership over them. Therefore
         * unique pointers are used.
         * @see Features
         */
        std::vector<std::unique_ptr<Features::ContinentalPlateModels::Temperature::Interface> > temperature_models;

        /**
         * A vector containing all the pointers to the composition models. This vector is
         * responsible for the features and has ownership over them. Therefore
         * unique pointers are used.
         * @see Features
         */
        std::vector<std::unique_ptr<Features::ContinentalPlateModels::Composition::Interface> > composition_models;

        /**
         * A vector containing all the pointers to the grains models. This vector is
         * responsible for the features and has ownership over them. Therefore
         * unique pointers are used.
         * @see Features
         */
        std::vector<std::unique_ptr<Features::ContinentalPlateModels::Grains::Interface> > grains_models;


        double min_depth;
        Objects::Surface min_depth_surface;
        double max_depth;
        Objects::Surface max_depth_surface;

    };


  } // namespace Features
} // namespace WorldBuilder

#endif
