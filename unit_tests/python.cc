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
