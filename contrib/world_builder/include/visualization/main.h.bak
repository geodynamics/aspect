/*
  Copyright (C) 2020 by the authors of the World Builder code.

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

#ifndef WORLD_BUILDER_VISUALIZATION_MAIN_H_
#define WORLD_BUILDER_VISUALIZATION_MAIN_H_

#include <vector>
#include <string>

void project_on_sphere(double, double &, double &, double &);

void lay_points(double x1, double y1, double z1,
                double x2, double y2, double z2,
                double x3, double y3, double z3,
                double x4, double y4, double z4,
                std::vector<double> &x, std::vector<double> &y, std::vector<double> &z,
                std::vector<bool> &hull, size_t level);

std::vector<std::string> get_command_line_options_vector(int argc, char **argv);

bool find_command_line_option(char **begin, char **end, const std::string &option);

#endif