/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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

#include <aspect/vofinterface/QCTrap.h>
#include <deal.II/base/geometry_info.h>

DEAL_II_NAMESPACE_OPEN

template <>
QCTrap<0>::QCTrap (const unsigned int)
  : Quadrature<0> (1)
{
  this->weights[0] = 1.0;
}

template <>
QCTrap<1>::QCTrap (const unsigned int n)
  : Quadrature<1> (n + 1)
{
  unsigned int i;
  double b_weight = 1.0 / n;

  if (n == 0)
    return;

  this->quadrature_points[0] = Point<1> (0.0);
  this->weights[0] = 0.5 * b_weight;
  for (i = 1; i < n; ++i)
    {
      this->quadrature_points[i] = Point<1> (i * b_weight);
      this->weights[i] = b_weight;
    }
  this->quadrature_points[n] = Point<1> (1.0);
  this->weights[n] = 0.5 * b_weight;
}

template <int dim>
QCTrap<dim>::QCTrap (const unsigned int n)
  : Quadrature<dim> (QCTrap<dim - 1> (n), QCTrap<1> (n))
{
}

template class QCTrap<2> ;
template class QCTrap<3> ;

DEAL_II_NAMESPACE_CLOSE

