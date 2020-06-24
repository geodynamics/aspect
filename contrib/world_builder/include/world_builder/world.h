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

#ifndef _world_builder_world_h
#define _world_builder_world_h

#include <random>

#include <world_builder/parameters.h>
#include <world_builder/grains.h>



namespace WorldBuilder
{


  namespace Features
  {
    class Interface;
  }

  class World
  {
    public:
      /**
       * Constructor. This constructor requires the parameter filename. Other parameters
       * are optional.
       * \param filename  a string with the location of
       * the world builder file to initialize the world.
       * \param has_output_dir a bool indicating whether the world builder is allowed to
       * write out information to a directly.
       * \param output_dir a string with the location of the directory where the world builder
       * is allowed to write information to if it is allowed by the bool has_output_dir.
       * \param random_number_seed a double containing a seed for the random number generator.
       * The world builder uses a deterministic random number generator for some plugins. This
       * is a deterministic random number generator on prorpose because even though you might
       * want to use random numbers to initialize some fields, the result should be reproducable.
       * Note that when the world builder is used in for example MPI programs you should supply
       * the world builder created each MPI process a different seed. You can use the MPI RANK
       * for this (seed is seed + MPI_RANK). Because the generator is deterministic (known and
       * documented algorithm), we can test the results and they should be the same even for different
       * compilers and machines.
       */
      World(std::string filename, bool has_output_dir = false, std::string output_dir = "", unsigned long random_number_seed = 1);

      /**
       * Destructor
       */
      ~World();


      /**
       * Describe what the world builder file should look like
       */
      static void declare_entries(Parameters &prm);

      /**
       * read in the world builder file
       */
      void parse_entries(Parameters &parameters);

      /**
       * Returns the temperature based on a 2d Cartesian point, the depth in the
       * model at that point and the gravity norm at that point.
       */
      double temperature(const std::array<double, 2> &point, const double depth, const double gravity_norm) const;

      /**
       * Returns the temperature based on a 3d Cartesian point, the depth in the
       * model at that point and the gravity norm at that point.
       */
      double temperature(const std::array<double, 3> &point, const double depth, const double gravity_norm) const;

      /**
       * Returns the composition value based on a 2d Cartesian point, the depth in
       * the model at that point and the gravity norm at that point.
       */
      double composition(const std::array<double, 2> &point, const double depth, const unsigned int composition_number) const;

      /**
       * Returns the composition value based on a 3d Cartesian point, the depth in
       * the model at that point and the gravity norm at that point.
       */
      double composition(const std::array<double, 3> &point, const double depth, const unsigned int composition_number) const;

      /**
       * Returns the grain orientations and sizes based on a 2d Cartesian point, the depth in
       * the model at that point and the gravity norm at that point.
       */
      WorldBuilder::grains grains(const std::array<double, 2> &point,
                                  const double depth,
                                  const unsigned int composition_number,
                                  size_t number_of_grains) const;

      /**
       * Returns the grain orientations and sizes based on a 3d Cartesian point, the depth in
       * the model at that point and the gravity norm at that point.
       */
      WorldBuilder::grains grains(const std::array<double, 3> &point,
                                  const double depth,
                                  const unsigned int composition_number,
                                  size_t number_of_grains) const;


      /**
       * Return a reference to the mt19937 random number.
       * The seed is provided to the world builder at construction.
       */
      std::mt19937 &get_random_number_engine();

      /**
       * This is the parameter class, which stores all the values loaded in
       * from the parameter file or which are set directly.
       */
      Parameters parameters;

      /**
       * Todo
       */
      std::vector<Point<2> > cross_section;

      /**
       * Todo
       */
      Point<2> surface_coord_conversions;

      /**
       * Todo
       */
      double potential_mantle_temperature;

      /**
       * Todo
       */
      double surface_temperature;

      /**
       * Todo
       */
      bool force_surface_temperature;

      /**
       * Todo
       */
      double thermal_expansion_coefficient;

      /**
       * Todo
       */
      double specific_heat;

      /**
       * Todo
       */
      double thermal_diffusivity;

      /**
       * Todo
       */
      double maximum_distance_between_coordinates;

      /**
       * Todo
       */
      std::string interpolation;


    private:
      /**
       * The minimum dimension. If cross section data is provided, it is set
       * to 2, which means the 2d function of temperature and composition can
       * be used. Otherwise it is set to 3, which means that they can't be
       * used.
       */
      unsigned int dim;


      /**
       * random number generator engine
       */
      std::mt19937 random_number_engine;



  };
}

#endif
