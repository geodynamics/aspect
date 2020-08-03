#!/bin/perl
#
# This script removes references to the std_cxx11 namespace of deal.II and
# replaces them with references to std, which is possible because we now
# require C++11 compliant compilers.

while (<>)
{
  s/cxx1x/cxx11/g;
  s:deal.II/base/std_cxx11/tuple.h:tuple:;
  s:deal.II/base/std_cxx11/array.h:array:;
  s:deal.II/base/std_cxx11/shared_ptr.h:memory:;
  s:deal.II/base/std_cxx11/unique_ptr.h:memory:;
  s:deal.II/base/std_cxx11/bind.h:functional:;
  s:deal.II/base/std_cxx11/function.h:functional:;

  s/std_cxx11::_/std::placeholders::_/g;

  s/std_cxx11/std/g;
  print;
}
