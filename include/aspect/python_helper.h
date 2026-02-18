/*
  Copyright (C) 2025 - 2026 by the authors of the ASPECT code.

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


#ifndef _aspect_python_helper_h
#define _aspect_python_helper_h

#include <aspect/global.h>
#include <deal.II/base/array_view.h>


#ifdef ASPECT_WITH_PYTHON

// Python does not like it if this macro is already defined. This happens at
// least in some versions of Trilinos and can trigger only with certain unity
// build options:
#ifdef HAVE_SYS_TIME_H
#  undef HAVE_SYS_TIME_H
#endif
#define PY_SSIZE_T_CLEAN
#include <Python.h>

// Declare that we are compatible with numpy 1.7 API and don't want any
// compile warnings for deprecated API calls:
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

// Define a unique symbol name for the numpy Python API for used when building
// the ASPECT executable. This will be used by numpy to declare a static member
// variable that will be initialized using the import_array() function (see below).
#ifndef PY_ARRAY_UNIQUE_SYMBOL
#  define PY_ARRAY_UNIQUE_SYMBOL ASPECT_ARRAY_API
#endif

// NumPY API requires a call to import_array() to be made before usage
// but only once per executable. We will define
// ASPECT_NUMPY_DEFINE_API in exactly one translation unit per
// executable/shared library (in main.cc and any plugin that uses
// Python). All other translation units including this header will
// therefore define NO_IMPORT_ARRAY that asks numpy to skip defining
// import_array().
#ifndef ASPECT_NUMPY_DEFINE_API
#  define NO_IMPORT_ARRAY
#endif

#include <numpy/arrayobject.h>

// Clean up any of the macros defined above. This is important if we
// are doing a unity build:
#ifdef ASPECT_NUMPY_DEFINE_API
#  undef ASPECT_NUMPY_DEFINE_API
#endif
#ifdef NO_IMPORT_ARRAY
#  undef NO_IMPORT_ARRAY
#endif

namespace aspect
{
  namespace PythonHelper
  {
    /**
     * Convert a std::vector<double> to a numpy array.
     *
     * The caller is responsible for decrementing the reference count of the returned PyObject once
     * it is no longer needed.
     * @param[in] vec The vector to convert.
     * @return A PyObject* to the numpy array inside a std::unique_ptr.
     */
    inline
    std::unique_ptr<PyObject, void(*)(PyObject *)> vector_to_numpy_object(const std::vector<double> &vec)
    {
      const npy_intp size = static_cast<npy_intp>(vec.size());

      auto deleter = [](PyObject *obj)
      {
        Py_DECREF(obj);
      };
      if (size == 0)
        return std::unique_ptr<PyObject, void(*)(PyObject *)>(PyArray_SimpleNew(1, &size, NPY_DOUBLE), deleter);
      double *data = const_cast<double *>(vec.data());
      return std::unique_ptr<PyObject, void(*)(PyObject *)>(PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, data), deleter);
    }



    /**
     * Access the contents of a numpy array (PyObject*) using an ArrayView<double>.
     * The PyObject must remain valid while the view is in use.
     */
    inline dealii::ArrayView<double> numpy_to_array_view(PyObject *obj)
    {
      AssertThrow(PyArray_Check(obj), dealii::ExcMessage("Expected a numpy array"));
      PyArrayObject *arr = reinterpret_cast<PyArrayObject *>(obj);
      return dealii::ArrayView<double>(static_cast<double *>(PyArray_DATA(arr)),
                                       static_cast<size_t>(PyArray_SIZE(arr)));
    }

  }
}
#endif


#endif
