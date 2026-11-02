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

#include <cstdlib>
#include <iostream>
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
  PyObject *pModule = PyImport_ImportModule("python_01");
  if (pModule == nullptr)
    {
      PyErr_Print();
      AssertThrow(false, dealii::ExcMessage("Failed to load Python module"));
    }

  std::vector<double> x = {1.0, 2.0, 3.0};
  auto arr = aspect::PythonHelper::vector_to_numpy_object(x);

  PyObject *pFunc = PyObject_GetAttrString(pModule, "update");
  Assert(pFunc != nullptr, dealii::ExcMessage("update() function not found."));
  PyObject *pArgs = PyTuple_Pack(2, PyFloat_FromDouble(42.0), arr.get());

  PyObject *result = PyObject_CallObject(pFunc, pArgs);
  Assert(result != nullptr, dealii::ExcMessage("Python call to update() failed."));
  double result_double = PyFloat_AsDouble(result);
  Assert(result_double == 21.0, dealii::ExcMessage("Incorrect result."));
  std::cout << "result: " << result_double << std::endl;

  Py_DECREF(pArgs);
  Py_DECREF(pFunc);

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
