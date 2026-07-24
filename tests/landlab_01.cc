/*
  Copyright (C) 2025 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/utilities.h>

#include <iostream>
#include <cfenv>

#define ASPECT_NUMPY_DEFINE_API
#include <aspect/python_helper.h>

// create a function that is run upon loading the plugin
int f()
{
  if (_import_array() < 0)
    {
      PyErr_Print();
      AssertThrow(false, dealii::ExcMessage("Numpy init failed!"));
    }

  PyRun_SimpleString("import sys; sys.path.append(\"" ASPECT_SOURCE_DIR "/tests\")");

  // avoid floating point exceptions in Landlab Python code:
#ifdef ASPECT_USE_FP_EXCEPTIONS
  fedisableexcept(FE_DIVBYZERO|FE_INVALID);
#endif
  PyObject *pModule = PyImport_ImportModule("landlab_01");
  if (pModule == nullptr)
    {
      PyErr_Print();
      AssertThrow(false, dealii::ExcMessage("Failed to load Python module"));
    }

  return 42;
}

template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &/*signals*/)
{
  f();
  std::cout << "exiting..." << std::endl;
  std::exit(0); // let's exit and not actually run ASPECT...
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
