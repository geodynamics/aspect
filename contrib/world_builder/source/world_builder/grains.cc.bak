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

#include "world_builder/grains.h"

namespace WorldBuilder
{

  grains::grains()


    = default;


  grains::grains(const std::vector<double> &vector,
                 const size_t number_of_grains,
                 const size_t start_entry)
  {
    sizes.resize(number_of_grains);
    for (unsigned int i_grain = 0; i_grain < number_of_grains; i_grain++)
      {
        sizes[i_grain] = vector[start_entry+i_grain];
      }

    rotation_matrices.resize(number_of_grains);
    for (unsigned int i_grain = 0; i_grain < number_of_grains; i_grain++)
      {
        rotation_matrices[i_grain][0][0] = vector[start_entry+number_of_grains+i_grain*9];
        rotation_matrices[i_grain][0][1] = vector[start_entry+number_of_grains+i_grain*9+1];
        rotation_matrices[i_grain][0][2] = vector[start_entry+number_of_grains+i_grain*9+2];
        rotation_matrices[i_grain][1][0] = vector[start_entry+number_of_grains+i_grain*9+3];
        rotation_matrices[i_grain][1][1] = vector[start_entry+number_of_grains+i_grain*9+4];
        rotation_matrices[i_grain][1][2] = vector[start_entry+number_of_grains+i_grain*9+5];
        rotation_matrices[i_grain][2][0] = vector[start_entry+number_of_grains+i_grain*9+6];
        rotation_matrices[i_grain][2][1] = vector[start_entry+number_of_grains+i_grain*9+7];
        rotation_matrices[i_grain][2][2] = vector[start_entry+number_of_grains+i_grain*9+8];
      }
  }


  void
  grains::unroll_into(std::vector<double> &vector, const size_t start_entry) const
  {
    const size_t number_of_grains = sizes.size();

    for (unsigned int i_grain = 0; i_grain < number_of_grains; i_grain++)
      {
        vector[start_entry+i_grain] = sizes[i_grain];
      }

    for (unsigned int i_grain = 0; i_grain < number_of_grains; i_grain++)
      {
        vector[start_entry+number_of_grains+i_grain*9]   = rotation_matrices[i_grain][0][0];
        vector[start_entry+number_of_grains+i_grain*9+1] = rotation_matrices[i_grain][0][1];
        vector[start_entry+number_of_grains+i_grain*9+2] = rotation_matrices[i_grain][0][2];
        vector[start_entry+number_of_grains+i_grain*9+3] = rotation_matrices[i_grain][1][0];
        vector[start_entry+number_of_grains+i_grain*9+4] = rotation_matrices[i_grain][1][1];
        vector[start_entry+number_of_grains+i_grain*9+5] = rotation_matrices[i_grain][1][2];
        vector[start_entry+number_of_grains+i_grain*9+6] = rotation_matrices[i_grain][2][0];
        vector[start_entry+number_of_grains+i_grain*9+7] = rotation_matrices[i_grain][2][1];
        vector[start_entry+number_of_grains+i_grain*9+8] = rotation_matrices[i_grain][2][2];
      }
  }
} // namespace WorldBuilder