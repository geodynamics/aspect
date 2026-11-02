// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#include <aspect/utilities.h>
#include <aspect/simulator_signals.h>

#include <iostream>
#include <cstdlib>
#include <cfenv>
#include <string>

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
  std::cout << "ok" << std::endl;
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
