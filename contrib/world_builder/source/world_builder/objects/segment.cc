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
#include <utility>

#include "world_builder/objects/segment.h"


namespace WorldBuilder
{

  namespace Features
  {
    namespace FaultModels
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
    }  // namespace FaultModels
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
  }  // namespace Features


  namespace Objects
  {
    // todo update function
    template<class A, class B, class C>
    Segment<A,B,C>::Segment(const double default_length_,
                            const WorldBuilder::Point<2> &default_thickness_,
                            const WorldBuilder::Point<2> &default_top_truncation_,
                            const WorldBuilder::Point<2> &default_angle_,
                            const std::vector<std::shared_ptr<A> > temperature_systems_,
                            const std::vector<std::shared_ptr<B> > composition_systems_,
                            const std::vector<std::shared_ptr<C> > grains_systems_)
      :
      value_length(default_length_),
      default_length(default_length_),
      value_thickness(default_thickness_),
      value_top_truncation(default_top_truncation_),
      value_angle(default_angle_),
      temperature_systems(std::move(temperature_systems_)),
      composition_systems(std::move(composition_systems_)),
      grains_systems(std::move(grains_systems_))
    {
    }

    template<class A, class B, class C>
    Segment<A,B,C>::Segment(Segment const &other)
      :
      value_length(other.value_length),
      default_length(other.default_length),
      value_thickness(other.value_thickness),
      value_top_truncation(other.value_top_truncation),
      value_angle(other.value_angle),
      temperature_systems(other.temperature_systems),
      composition_systems(other.composition_systems),
      grains_systems(other.grains_systems)
    {
    }

    template<class A, class B, class C>
    Segment<A,B,C>::~Segment ()
      = default;


    /**
    * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
    * Note that the variable with this name has to be loaded before this function is called.
    */
    template class
    Segment<Features::SubductingPlateModels::Temperature::Interface,Features::SubductingPlateModels::Composition::Interface,Features::SubductingPlateModels::Grains::Interface>;

    /**
    * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
    * Note that the variable with this name has to be loaded before this function is called.
    */
    template class
    Segment<Features::FaultModels::Temperature::Interface,Features::FaultModels::Composition::Interface,Features::FaultModels::Grains::Interface>;

    /**
    * Todo: Returns a vector of pointers to the Point<3> Type based on the provided name.
    * Note that the variable with this name has to be loaded before this function is called.
    */
    //template class
    //Segment<char,char>;

  } // namespace Objects
} // namespace WorldBuilder