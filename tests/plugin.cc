// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2022 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

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
