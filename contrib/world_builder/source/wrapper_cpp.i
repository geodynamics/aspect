%begin %{
#ifdef _MSC_VER
#define SWIG_PYTHON_INTERPRETER_NO_DEBUG
#endif
%}
%include <stl.i>
%module gwb
%{
#include "../include/world_builder/wrapper_cpp.h"
%}
%include "../include/world_builder/wrapper_cpp.h"
