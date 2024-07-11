/*
  Copyright (C) 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/utilities.h>

#include <iostream>

void print_string_vector(std::vector<std::string> data)
{
  for (std::vector<std::string>::iterator it = data.begin(); it != data.end(); ++it)
    std::cout << "\t" << *it << std::endl;
}

int f()
{
  using namespace aspect::Utilities;

  std::vector<std::string> base_declarations;
  base_declarations.push_back("vector(  vec  )");
  base_declarations.push_back("tensor(  ten )");

  std::vector<std::string> var_decl_2 = expand_dimensional_variable_names<2> (base_declarations);
  std::vector<std::string> var_decl_3 = expand_dimensional_variable_names<3> (base_declarations);

  std::cout << "Dimension 2:" << std::endl;
  print_string_vector(var_decl_2);

  std::cout << "Dimension 3:" << std::endl;
  print_string_vector(var_decl_3);

  exit(0);

  // dummy return for compilation warnings

  return 0;
}

int ignored_variable = f();
