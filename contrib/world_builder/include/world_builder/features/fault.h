/*
  Copyright (C) 2018 - 2020 by the authors of the World Builder code.

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

#ifndef _world_feature_features_fault_h
#define _world_feature_features_fault_h

#include <world_builder/features/interface.h>
#include <world_builder/world.h>
#include <world_builder/types/segment.h>

#include <world_builder/features/fault_models/temperature/interface.h>
#include <world_builder/features/fault_models/composition/interface.h>
#include <world_builder/features/fault_models/grains/interface.h>


namespace WorldBuilder
{
  namespace Features
  {
    /**
     * This class represents a fault and can implement submodules
     * for temperature and composition. These submodules determine what
     * the returned temperature or composition of the temperature and composition
     * functions of this class will be.
     */
    class Fault : public Interface
    {
      public:
        /**
         * constructor
         */
        Fault(WorldBuilder::World *world);

        /**
         * Destructor
         */
        ~Fault();

        /**
         * declare and read in the world builder file into the parameters class
         */
        static
        void declare_entries(Parameters &prm,
                             const std::string &parent_name = "",
                             const std::vector<std::string> &required_entries = {});

        /**
         * declare and read in the world builder file into the parameters class
         */
        void parse_entries(Parameters &prm) override final;

        /**
         * Returns a temperature based on the given position, depth in the model,
         * gravity and current temperature.
         */
        double temperature(const Point<3> &position,
                           const double depth,
                           const double gravity,
                           double temperature) const override final;

        /**
         * Returns a value for the requests composition (0 is not present,
         * 1 is present) based on the given position, depth in the model,
         * the composition which is being requested and the current value
         * of that composition at this location and depth.
         */
        double composition(const Point<3> &position,
                           const double depth,
                           const unsigned int composition_number,
                           double composition_value) const override final;

        /**
         * Takes the list of the fault coordinates and computes a bounding box 
         * based on the minimum and maximum x,y values in the list, and an
         * additional buffer zone computed using the maximum fault length.
         * Currenly, only implemented for the spherical coordinates.
         */

        std::vector<Point<2> >
        get_bounding_polygon (const std::vector<Point<2> > &point_list) const;

        /**
         * Takes the list of coordinates defining the bounding box, and returns 
         * true if the point lies inside the box and false otherwise.
         */
        bool bounding_polygon_contains_point (const std::vector<Point<2> > &point_list,
                                              const Point<2> &point) const;

        /**
         * Returns a grains (rotation matrix and grain size)
         * based on the given position, depth in the model,
         * the composition which is being requested and the current value
         * of that composition at this location and depth.
         */
        virtual
        WorldBuilder::grains
        grains(const Point<3> &position,
               const double depth,
               const unsigned int composition_number,
               WorldBuilder::grains value) const override final;



      private:
        std::vector<std::shared_ptr<Features::FaultModels::Temperature::Interface> > default_temperature_models;
        std::vector<std::shared_ptr<Features::FaultModels::Composition::Interface>  > default_composition_models;
        std::vector<std::shared_ptr<Features::FaultModels::Grains::Interface>  > default_grains_models;

        std::vector<Objects::Segment<Features::FaultModels::Temperature::Interface,
            Features::FaultModels::Composition::Interface,
            Features::FaultModels::Grains::Interface> > default_segment_vector;

        std::vector< std::vector<Objects::Segment<Features::FaultModels::Temperature::Interface,
            Features::FaultModels::Composition::Interface,
            Features::FaultModels::Grains::Interface> > > sections_segment_vector;

        // This vector stores segments to this coordiante/section.
        //First used (raw) pointers to the segment relevant to this coordinate/section,
        // but I do not trust it won't fail when memory is moved. So storing the all the data now.
        std::vector<std::vector<Objects::Segment<Features::FaultModels::Temperature::Interface,
            Features::FaultModels::Composition::Interface,
            Features::FaultModels::Grains::Interface> > > segment_vector;

        // todo: the memory of this can be greatly improved by
        // or using a plugin system for the submodules, or
        // putting the variables in a union. Although the memory
        // used by this program is in real cases expected to be
        // insignificant compared to what a calling program may
        // use, a smaller amount of memory used in here, could
        // theoretically speed up the computation, because more
        // relevant data could be stored in the cache. But this
        // not a urgent problem, and would require testing.

        /**
         * This variable stores the depth at which the fault
         * starts. It makes this depth effectively the surface
         * of the model for the fault.
         */
        double starting_depth;

        /**
         * The depth which below the fault may no longer
         * be present. This can not only help setting up models with
         * less effort, but can also improve performance, because
         * the algorithm doesn't have to search in locations below
         * this depth.
         */
        double maximum_depth;

        /**
         * A point on the surface to which the fault dips.
         */
        WorldBuilder::Point<2> reference_point;

        std::vector<std::vector<double> > slab_segment_lengths;
        std::vector<std::vector<Point<2> > > slab_segment_thickness;
        std::vector<std::vector<Point<2> > > slab_segment_top_truncation;
        std::vector<std::vector<Point<2> > > slab_segment_angles;
        std::vector<double> total_slab_length;
        double maximum_total_slab_length;
        double maximum_slab_thickness;


    };
  }
}

#endif
