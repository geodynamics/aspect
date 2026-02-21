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

#ifndef WORLD_BUILDER_WRAPPER_C_H
#define WORLD_BUILDER_WRAPPER_C_H

#include <stdbool.h>

/** @file
 * This file contains the C API of the WorldBuilder library.
 *
 * The C API is used to make the WorldBuilder library available to programs written
 * in C. The C API is a set of functions that can be called from C code. The
 * functions in the C API call the corresponding C++ functions in the WorldBuilder
 * library. Check the content of this file for a list of the available functions and
 * their documentation.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * This function creates an object of the world builder and returns a pointer
 * to it. This pointer can then be used to call the temperature and composition
 * functions. When done call the release world function to destroy the object.
 */
void create_world(void **ptr_ptr_world, const char *world_builder_file, const bool *has_output_dir, const char *output_dir, const unsigned long random_number_seed);



/**
 * Return the size of the output vector returned by the properties function for a given properties vector.
 *
 * @param properties The properties parameter from the properties function. See the documentation of that
 * function for more info.
 * @param n_properties number of properties
 * @return unsigned int Return the size of the output vector returned by the properties function for a given properties vector.
 */
unsigned int properties_output_size(void *ptr_ptr_world,
                                    const unsigned int properties[][3],
                                    const unsigned int n_properties);

/**
 * @brief This function returns 2D properties.
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
 *
 * @param ptr_ptr_world a pointer to the world
 * @param x The x position of the point
 * @param z The z position of the point
 * @param depth The depth of the point
 * @param properties an array of properties, which each property is an array of three integers.
 * @param n_properties number of properties.
 * @param values are the return values as a pointer to an array of doubles
 */
void properties_2d(void *ptr_ptr_world,
                   const double x,
                   const double z,
                   const double depth,
                   const unsigned int properties[][3],
                   const unsigned int n_properties,
                   double values[]);


/**
 * @brief This function returns 3D properties.
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
 *
 * @param ptr_ptr_world a pointer to the world
 * @param x The x position of the point
 * @param x The y position of the point
 * @param z The z position of the point
 * @param depth The depth of the point
 * @param properties an array of properties, which each property is an array of three integers.
 * @param n_properties number of properties.
 * @param values are the return values as a pointer to an array of doubles
 */
void properties_3d(void *ptr_ptr_world,
                   const double x,
                   const double y,
                   const double z,
                   const double depth,
                   const unsigned int properties[][3],
                   const unsigned int n_properties,
                   double values[]);

/**
 * This function return the temperature at a specific location given x, z, depth and
 * gravity.
 * Note: gravity value is no longer used, instead use the gravity model from the input file.
 */
void temperature_2d(void *ptr_ptr_world, double x, double z, double depth, double *temperature);

/**
 * This function return the temperature at a specific location given x, y, z, depth and
 * gravity.
 */
void temperature_3d(void *ptr_ptr_world, double x, double y, double z, double depth, double *temperature);

/**
 * This function return the composition at a specific location given x, z, depth and
 * composition number.
 */
void composition_2d(void *ptr_ptr_world, double x, double z, double depth, unsigned int composition_number, double *composition);

/**
 * This function return the composition at a specific location given x, y, z, depth and
 * composition number.
 */
void composition_3d(void *ptr_ptr_world, double x, double y, double z, double depth, unsigned int composition_number, double *composition);

/**
 * The destructor for the world builder class. Call this function when done with the
 * world builder.
 */
void release_world(void *ptr_ptr_world);

#ifdef __cplusplus
}
#endif

#endif
