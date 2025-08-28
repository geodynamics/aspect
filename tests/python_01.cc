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

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

PyObject *vector_to_numpy(const std::vector<double> &vec)
{
  npy_intp size = vec.size();
  double *data = const_cast<double *>(vec.data());
  return PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data);
}

// create a function that is run upon loading the plugin
int f()
{
  // required for Numpy interop in each translation unit
  if (_import_array() < 0)
    AssertThrow(false, dealii::ExcMessage("Numpy init failed!"));

  PyRun_SimpleString("import sys; sys.path.append(\"" ASPECT_SOURCE_DIR "/tests\")");
  PyObject *pModule = PyImport_ImportModule("python_01");
  AssertThrow(pModule != nullptr, dealii::ExcMessage("Failed to load Python module"));

  std::vector<double> x = {1.0, 2.0, 3.0};
  PyObject *arr = vector_to_numpy(x);

  PyObject *pFunc = PyObject_GetAttrString(pModule, "update");
  Assert(pFunc != nullptr, dealii::ExcMessage("update() function not found."));
  PyObject *pArgs = PyTuple_Pack(2, PyFloat_FromDouble(42.0), arr);

  PyObject *result = PyObject_CallObject(pFunc, pArgs);
  Assert(result != nullptr, dealii::ExcMessage("Python call to update() failed."));
  double result_double = PyFloat_AsDouble(result);
  Assert(result_double == 21.0, dealii::ExcMessage("Incorrect result."));
  std::cout << "result: " << result_double << std::endl;


  Py_DECREF(arr);
  Py_DECREF(pArgs);
  Py_DECREF(pFunc);

  return 42;
}

template <int dim>
void signal_connector (aspect::SimulatorSignals<dim> &/*signals*/)
{
  f();
  std::cout << "exiting..." << std::endl;
  std::exit(0); // let's exist and not actually run ASPECT...
}

ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                  signal_connector<3>)
