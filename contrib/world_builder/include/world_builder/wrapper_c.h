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
