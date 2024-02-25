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

#ifndef WORLD_BUILDER_WORLD_H
#define WORLD_BUILDER_WORLD_H

#include "world_builder/grains.h"
#include "world_builder/parameters.h"
#include "world_builder/utilities.h"
#include "world_builder/objects/distance_from_surface.h"

#include <random>

/**
* The global namespace for the Geodynamic World Builder
*/
namespace WorldBuilder
{

  namespace Features
  {
    class Interface;
  } // namespace Features

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
       * want to use random numbers to initialize some fields, the result should be reproducible.
       * Note that when the world builder is used in for example MPI programs you should supply
       * the world builder created each MPI process a different seed. You can use the MPI RANK
       * for this (seed is seed + MPI_RANK). Because the generator is deterministic (known and
       * documented algorithm), we can test the results and they should be the same even for different
       * compilers and machines.
       */
      World(std::string filename, bool has_output_dir = false, const std::string &output_dir = "", unsigned long random_number_seed = 1, const bool limit_debug_consistency_checks = true);

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
      void parse_entries(Parameters &prm);

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
       */
      std::vector<double> properties(const std::array<double, 2> &point,
                                     const double depth,
                                     const std::vector<std::array<unsigned int,3>> &properties) const;

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
       * number is not used. So a composition query asking about composition 1 looks like this:
       * {2,1,0}. A composition query prodoces one entry in the output vector.
       *
       * Grains are identified by 3. The second entry is the grain composition number and the third
       * entry is the number of grains. A query about the grains, where it asks about composition 1
       * (for example enstatite) and 500 grains, looks like this: {3,1,500}.
       * A composition query prodoces n_grains*10 entries in the output vector. The first n_grains
       * entries are the sizes of all the grains, and the other 9 entries are sets of rotation
       * matrices. The rotation matrix entries are ordered [0][0],[0][1],[0][2],[1][0],[1][1],etc.
       *
       * The tag is identified by 4 and no extra information is needed. So the tag
       * input usually looks like {4,0,0}. A tag query prodoces one entry in the output
       * vector, representing the index of the tag of the last/dominant feature.
       */
      std::vector<double> properties(const std::array<double, 3> &point,
                                     const double depth,
                                     const std::vector<std::array<unsigned int,3>> &properties) const;

      /**
       * Returns the temperature based on a 2d Cartesian point, the depth in the
       * model at that point and the gravity norm at that point.
       */
      double temperature(const std::array<double, 2> &point, const double depth) const;

      /**
       * Returns the temperature based on a 3d Cartesian point, the depth in the
       * model at that point and the gravity norm at that point.
       */
      double temperature(const std::array<double, 3> &point, const double depth) const;

      /**
       * Returns the temperature based on a 2d Cartesian point, the depth in the
       * model at that point and the gravity norm at that point.
       * Note: gravity norm is no longer used, instead use the gravity model from the input file.
       */
      [[deprecated("Replaced by a temperature function without the gravity. This function will be removed in future versions.")]]
      double temperature(const std::array<double, 2> &point, const double depth, const double gravity_norm) const;

      /**
       * Returns the temperature based on a 3d Cartesian point, the depth in the
       * model at that point and the gravity norm at that point.
       * Note: gravity norm is no longer used, instead use the gravity model from the input file.
       */
      [[deprecated("Replaced by a temperature function without the gravity. This function will be removed in future versions.")]]
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
       * Returns a PlaneDistances object that has the distance from and along a feature plane,
       * calculated from the coordinates and the depth of the point.
       \param point the coordinates in the cartesian geometry
       \param depth the depth of the point
       \param name the name of the feature (i.e. the string provided to the key word "name" in the wb file)
      */
      Objects::PlaneDistances
      distance_to_plane(const std::array<double, 3> &point,
                        const double depth,
                        const std::string &name) const;

      /**
       * The MPI rank. Set to zero if MPI is not available.
       */
      int MPI_RANK;

      /**
       * The MPI size. Set to one if MPI is not available.
       */
      int MPI_SIZE;

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

      /**
       * A list of all the feature tags.
       */
      std::vector<std::string> feature_tags;

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

      /**
       * limits some of the consistency checks in debug mode.
       * Current only prevents a check whether depth in spherical
       * coordinates is consistent with the computed depth from
       * x,y,z and provided radius.
       * Note: Recommended to keep it at false, unless you know what you are doing.
       */
      bool limit_debug_consistency_checks;



  };
} // namespace WorldBuilder

#endif
