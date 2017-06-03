/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/
#include <aspect/global.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/point.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/signaling_nan.h>

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>

#include <fstream>
#include <string>
#include <locale>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Plugins
  {
    InterfaceBase::~InterfaceBase ()
    {}



    void
    InterfaceBase::initialize ()
    {}



    void
    InterfaceBase::update ()
    {}



    void
    InterfaceBase::
    declare_parameters (dealii::ParameterHandler &)
    {}



    void
    InterfaceBase::parse_parameters (dealii::ParameterHandler &)
    {}

  }
}
