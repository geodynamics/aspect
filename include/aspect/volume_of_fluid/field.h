// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_volume_of_fluid_field_h
#define _aspect_volume_of_fluid_field_h

#include <aspect/fe_variable_collection.h>

namespace aspect
{
  namespace VolumeOfFluid
  {
    /**
     * A structure that contains enum values that identify type of input data
     * to allow robust use of sub-mesh scale input that remain valid over
     * multiple mesh sizes.
     */
    struct VolumeOfFluidInputType
    {
      enum Kind
      {
        /**
         * Input data is a value between 0 and 1 at all points.
         */
        composition,
        /**
         * Input data is an interface defined by a signed distance level set
         * with positive value indicating fluid presence. IE the function has
         * gradient 1 almost everywhere, is positive where the fluid is, and is
         * zero on the fluid interface
         */
        level_set
      };
    };
  }

  /**
   * Structure to package the relevant data (both state and cached) in a single
   * location for access.
   */
  template <int dim>
  struct VolumeOfFluidField
  {
    /**
     * Initialize the structure with FEVariables to hold the required
     * information that must be available on all cells.
     */
    VolumeOfFluidField(const FEVariable<dim> &volume_fraction,
                       const FEVariable<dim> &reconstruction,
                       const FEVariable<dim> &level_set,
                       const unsigned int composition_index);

    /**
     * Field to hold the current volume fraction.
     */
    const FEVariable<dim> &volume_fraction;

    /**
     * Field to hold the cached interface reconstruction.
     */
    const FEVariable<dim> &reconstruction;

    /**
     * Field to expose reconstructed interface as a zero-contour to output in
     * visualization plugin.
     */
    const FEVariable<dim> &level_set;

    /**
     * Field index of the associated composition field
     */
    const unsigned int composition_index;

  };
}

#endif
