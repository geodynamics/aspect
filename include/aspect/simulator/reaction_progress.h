// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#ifndef _aspect_simulator_reaction_progress_h
#define _aspect_simulator_reaction_progress_h

#include <algorithm>

namespace aspect
{
  namespace internal
  {
    /**
     * Apply a fixed-step reaction-progress increment, bounded to [0, 1],
     * and return the increment actually applied.
     */
    inline double
    apply_reaction_progress_increment(double &value,
                                      const double increment)
    {
      const double old_value = value;
      value = std::clamp(old_value + increment, 0.0, 1.0);
      return value - old_value;
    }
  }
}

#endif
