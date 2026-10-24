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

#include "common.h"
#include <aspect/utilities.h>

#include <aspect/python_helper.h>

#ifdef ASPECT_WITH_PYTHON

TEST_CASE("Utilities::python-string")
{
  using namespace dealii;

  PyObject *code_obj = Py_CompileString("1+1", "<string>", Py_eval_input);
  PyObject *main_module = PyImport_AddModule("__main__");
  PyObject *global_dict = PyModule_GetDict(main_module);
  PyObject *local_dict = PyDict_New();
  PyObject *obj = PyEval_EvalCode(code_obj, global_dict, local_dict);
  long result = PyLong_AsLong(obj);
  REQUIRE(result == 2);
}

TEST_CASE("Utilities::python-numpy-array1")
{
  using namespace dealii;

  const std::vector<double> x = {1.0, 2.0, 3.0};
  auto arr = aspect::PythonHelper::vector_to_numpy_object(x);
  ArrayView<double> view = aspect::PythonHelper::numpy_to_array_view(arr.get());
  REQUIRE(view.size() == 3);
  REQUIRE(view[0] == 1.0);
  REQUIRE(view[1] == 2.0);
  REQUIRE(view[2] == 3.0);
}

#endif
