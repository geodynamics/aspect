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

// make sure we can include deal.II and aspect files
#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <deal.II/grid/tria.h>

#include <iostream>
#include <typeinfo>

// create a function that is run upon loading the plugin
// and that produces some output
int f()
{
  std::cout << aspect::Utilities::expand_ASPECT_SOURCE_DIR("srcdir='$ASPECT_SOURCE_DIR'") << std::endl;

  std::cout << typeid(dealii::Triangulation<2>).name() << ' '
            << typeid(aspect::Simulator<2>).name() << std::endl;
  return 42;
}

// run this function by initializing a global variable by it
int i = f();
