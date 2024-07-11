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

#ifndef WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_H
#define WORLD_BUILDER_FEATURES_SUBDUCTING_PLATE_H


#include "world_builder/features/subducting_plate_models/composition/interface.h"
#include "world_builder/features/subducting_plate_models/grains/interface.h"
#include "world_builder/features/subducting_plate_models/temperature/interface.h"
#include "world_builder/objects/segment.h"
#include "world_builder/bounding_box.h"
#include "world_builder/objects/distance_from_surface.h"


namespace WorldBuilder
{
  class Parameters;
  class World;

  namespace Features
  {
    using namespace FeatureUtilities;
    namespace SubductingPlateModels
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
    }  // namespace SubductingPlateModels

    /**
     * This class represents a subducting plate and can implement submodules
     * for temperature and composition. These submodules determine what
     * the returned temperature or composition of the temperature and composition
     * functions of this class will be.
     */
    class SubductingPlate final: public Interface
    {
      public:
        /**
         * constructor
         */
        SubductingPlate(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~SubductingPlate() override final;

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
         * Returns the bounding points for a BoundingBox object using two extreme points in all the surface
         * coordinates and an additional buffer zone that accounts for the fault thickness and length. The first and second
         * points correspond to the lower left and the upper right corners of the bounding box, respectively (see the
         * documentation in include/bounding_box.h).
         * For the spherical system, the buffer zone along the longitudal direction is calculated using the
         * corresponding latitude points.
         */
        const BoundingBox<2>  &get_surface_bounding_box () const;

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

        /**
        * Returns a PlaneDistances object that has the distance from and along a subducting plate plane,
        * calculated from the coordinates and the depth of the point.
        */
        Objects::PlaneDistances
        distance_to_feature_plane(const Point<3> &position_in_cartesian_coordinates,
                                  const Objects::NaturalCoordinate &position_in_natural_coordinates,
                                  const double depth) const override;


      private:
        std::vector<std::shared_ptr<Features::SubductingPlateModels::Temperature::Interface> > default_temperature_models;
        std::vector<std::shared_ptr<Features::SubductingPlateModels::Composition::Interface>  > default_composition_models;
        std::vector<std::shared_ptr<Features::SubductingPlateModels::Grains::Interface>  > default_grains_models;

        std::vector<Objects::Segment<Features::SubductingPlateModels::Temperature::Interface,
            Features::SubductingPlateModels::Composition::Interface,
            Features::SubductingPlateModels::Grains::Interface> > default_segment_vector;

        std::vector< std::vector<Objects::Segment<Features::SubductingPlateModels::Temperature::Interface,
            Features::SubductingPlateModels::Composition::Interface,
            Features::SubductingPlateModels::Grains::Interface> > > sections_segment_vector;

        // This vector stores segments to this coordinate/section.
        //First used (raw) pointers to the segment relevant to this coordinate/section,
        // but I do not trust it won't fail when memory is moved. So storing the all the data now.
        std::vector<std::vector<Objects::Segment<Features::SubductingPlateModels::Temperature::Interface,
            Features::SubductingPlateModels::Composition::Interface,
            Features::SubductingPlateModels::Grains::Interface> > > segment_vector;

        // todo: the memory of this can be greatly improved by
        // or using a plugin system for the submodules, or
        // putting the variables in a union. Although the memory
        // used by this program is in real cases expected to be
        // Insignificant compared to what a calling program may
        // use, a smaller amount of memory used in here, could
        // theoretically speed up the computation, because more
        // relevant data could be stored in the cache. But this
        // not a urgent problem, and would require testing.

        /**
         * This variable stores the depth at which the subducting
         * plate starts. It makes this depth effectively the surface
         * of the model for the slab.
         */
        double starting_depth;

        /**
         * The depth which below the subducting plate may no longer
         * be present. This can not only help setting up models with
         * less effort, but can also improve performance, because
         * the algorithm doesn't have to search in locations below
         * this depth.
         */
        double maximum_depth;

        /**
         * Stores the bounding points for a BoundingBox object using two extreme points in all the surface
         * coordinates and an additional buffer zone that accounts for the fault thickness and length. The first and second
         * points correspond to the lower left and the upper right corners of the bounding box, respectively (see the
         * documentation in include/bounding_box.h).
         * For the spherical system, the buffer zone along the longitudal direction is calculated using the
         * corresponding latitude points.
         */
        BoundingBox<2> surface_bounding_box;

        /**
         * A point on the surface to which the subducting plates subduct.
         */
        Point<2> reference_point;

        std::vector<std::vector<double> > slab_segment_lengths;
        std::vector<std::vector<Point<2> > > slab_segment_thickness;
        std::vector<std::vector<Point<2> > > slab_segment_top_truncation;
        std::vector<std::vector<Point<2> > > slab_segment_angles;
        std::vector<double> total_slab_length;
        double maximum_total_slab_length;
        double maximum_slab_thickness;

        double min_along_x;
        double max_along_x;
        double min_along_y;
        double max_along_y;
        double min_lat_cos_inv;
        double max_lat_cos_inv;
        double buffer_around_slab_cartesian;
    };
  } // namespace Features
} // namespace WorldBuilder

#endif
