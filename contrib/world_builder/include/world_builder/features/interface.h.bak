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

#ifndef WORLD_BUILDER_FEATURES_INTERFACE_H
#define WORLD_BUILDER_FEATURES_INTERFACE_H


#include "world_builder/grains.h"
#include "world_builder/utilities.h"
#include "world_builder/objects/distance_from_surface.h"

namespace WorldBuilder
{
  class World;
  class Parameters;


  namespace Features
  {
    class ObjectFactory;

    /**
     * This class is an interface for the specific plate tectonic feature classes,
     * such as continental plate, oceanic plate and subduction zone.
     */
    class Interface
    {
      public:
        /**
         * constructor
         */
        Interface();

        /**
         * Destructor
         */
        virtual
        ~Interface();

        /**
         * helper function to parse coordinates.
         */
        void
        get_coordinates(const std::string &name,
                        Parameters &prm,
                        const CoordinateSystem coordinate_system);

        /**
         * declare and read in the world builder file into the parameters class
         */
        static
        void declare_entries(Parameters &prm,
                             const std::string &parent_name,
                             const std::vector<std::string> &required_entries);

        /**
         * declare and read in the world builder file into the parameters class
         */
        virtual
        void parse_entries(Parameters &prm) = 0;


        /**
         * takes a set of properties and a position and return a new set of properties
         */
        virtual
        void properties(const Point<3> &position_in_cartesian_coordinates,
                        const Objects::NaturalCoordinate &position_in_natural_coordinates,
                        const double depth,
                        const std::vector<std::array<unsigned int,3>> &properties,
                        const double gravity,
                        const std::vector<size_t> &entry_in_output,
                        std::vector<double> &output) const = 0;

        /**
         * A function to register a new type. This is part of the automatic
         * registration of the object factory.
         */
        static void registerType(const std::string &name,
                                 void ( * /*declare_entries*/)(Parameters &, const std::string &, const std::vector<std::string> &required_entries),
                                 void ( *make_snippet)(Parameters &),
                                 ObjectFactory *factory);

        std::string get_name() const
        {
          return name;
        };


        /**
         * A function to create a new type. This is part of the automatic
         * registration of the object factory.
         */
        static std::unique_ptr<Interface> create(const std::string &name, WorldBuilder::World *world);

        /**
        * Returns a PlaneDistances object that has the distance from and along a feature plane,
        * calculated from the coordinates and the depth of the point.
        */
        virtual
        Objects::PlaneDistances
        distance_to_feature_plane(const Point<3> &position_in_cartesian_coordinates,
                                  const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                  const double depth) const;


      protected:
        /**
         * A pointer to the world class to retrieve variables.
         */
        WorldBuilder::World *world;

        /**
         * The name of the feature type.
         */
        std::string name;

        /**
         * The index of the tag for this feature.
         * This corresponds to the index in the feature_tags variable which is store in the World.
         */
        size_t tag_index;

        /**
         * The type of interpolation used to get the line position between the points.
         */
        WorldBuilder::Utilities::InterpolationType interpolation_type;

        /**
         * number of original coordinates, before adding
         * more automatically.
         */
        std::size_t original_number_of_coordinates;

        /**
         * The coordinates at the surface of the feature
         */
        std::vector<Point<2> > coordinates;

        /**
         * The x and y spline
         */
        WorldBuilder::Objects::BezierCurve bezier_curve;


        /**
         * The name of the temperature submodule used by this feature.
         */
        std::string temperature_submodule_name;

        /**
         * The name of the composition submodule used by this feature.
         */
        std::string composition_submodule_name;


      private:
        static std::map<std::string, ObjectFactory *> &get_factory_map()
        {
          static std::map<std::string, ObjectFactory *> factories;
          return factories;
        }

        static std::map<std::string, void ( *)(Parameters &,
                                               const std::string &,
                                               const std::vector<std::string>& required_entries)> &get_declare_map()
        {
          static std::map<std::string, void ( *)(Parameters &,
                                                 const std::string &,
                                                 const std::vector<std::string>& required_entries)> declares;
          return declares;
        }

        static std::map<std::string, void ( *)(Parameters &)> &get_snippet_map()
        {
          static std::map<std::string, void ( *)(Parameters &)> declares;
          return declares;
        }
    };


    /**
     * A class to create new objects
     */
    class ObjectFactory
    {
      public:
        virtual std::unique_ptr<Interface> create(World *world) = 0;
    };

    /**
     * A macro which should be in every derived cpp file to automatically
     * register it. Because this is a library, we need some extra measures
     * to ensure that the static variable is actually initialized.
     */
#define WB_REGISTER_FEATURE(klass,name) \
  class klass##Factory : public ObjectFactory { \
    public: \
      klass##Factory() \
      { \
        Interface::registerType(#name, klass::declare_entries, klass::make_snippet, this); \
      } \
      std::unique_ptr<Interface> create(World *world) override final { \
        return std::unique_ptr<Interface>(new klass(world)); \
      } \
  }; \
  static klass##Factory global_##klass##Factory;

  } // namespace Features
} // namespace WorldBuilder

#endif
