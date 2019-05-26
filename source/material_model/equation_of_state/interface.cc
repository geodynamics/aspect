/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    EquationOfStateOutputs<dim>::EquationOfStateOutputs(const unsigned int n_points)
      :
      densities(n_points, numbers::signaling_nan<double>()),
      thermal_expansion_coefficients(n_points, numbers::signaling_nan<double>()),
      specific_heat(n_points, numbers::signaling_nan<double>()),
      compressibilities(n_points, numbers::signaling_nan<double>()),
      entropy_derivative_pressure(n_points, numbers::signaling_nan<double>()),
      entropy_derivative_temperature(n_points, numbers::signaling_nan<double>())
    {}
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  template struct EquationOfStateOutputs<dim>;
    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
