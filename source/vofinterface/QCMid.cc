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

#include <aspect/vofinterface/QCMid.h>
#include <deal.II/base/geometry_info.h>

DEAL_II_NAMESPACE_OPEN

template <>
QCMid<0>::QCMid (const unsigned int)
  : Quadrature<0> (1)
{
  this->weights[0] = 1.0;
}

template <>
QCMid<1>::QCMid (const unsigned int n)
  : Quadrature<1> (n)
{
  unsigned int i;
  double h = 1.0 / n;

  if (n == 0)
    return;

  for (i = 0; i < n; ++i)
    {
      this->quadrature_points[i] = Point<1> ((i+0.5) * h);
      this->weights[i] = h;
    }
}

template <int dim>
QCMid<dim>::QCMid (const unsigned int n)
  : Quadrature<dim> (QCMid<dim - 1> (n), QCMid<1> (n))
{
}

template class QCMid<2> ;
template class QCMid<3> ;

DEAL_II_NAMESPACE_CLOSE

