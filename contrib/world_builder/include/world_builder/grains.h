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

#ifndef WORLD_BUILDER_GRAINS_H
#define WORLD_BUILDER_GRAINS_H

#include <array>
#include <ostream>
#include <vector>

namespace WorldBuilder
{
  /**
   * This is a simple structure to store information about grains.
   * The advantage of storing all grains in separate vectors, compared
   * to having a vector of individual grains, is that the vectors can
   * be empty if the information is not needed.
   */
  struct grains
  {
    grains();

    grains(const std::vector<double> &vector,
           const size_t number_of_grains,
           const size_t start_entry = 0);

    void unroll_into(std::vector<double> &vector,
                     const size_t start_entry = 0) const;

    // The sizes of the grains
    std::vector<double> sizes;

    // the rotation matrices of the latices of the grains.
    // todo: convention.
    std::vector<std::array<std::array<double,3>,3> > rotation_matrices;

    friend std::ostream &operator<<(std::ostream &os, const grains &grains)
    {
      for (unsigned int i = 0; i < grains.sizes.size(); ++i)
        {
          os << i << ": s=" << grains.sizes[i] << ", R="
             << grains.rotation_matrices[i][0][0] << " " << grains.rotation_matrices[i][0][1] << " " << grains.rotation_matrices[i][0][2] << " "
             << grains.rotation_matrices[i][1][0] << " " << grains.rotation_matrices[i][1][1] << " " << grains.rotation_matrices[i][1][2] << " "
             << grains.rotation_matrices[i][2][0] << " " << grains.rotation_matrices[i][2][1] << " " << grains.rotation_matrices[i][2][2] << " ";
        }
      return os;
    }
  };

} // namespace WorldBuilder

#endif
