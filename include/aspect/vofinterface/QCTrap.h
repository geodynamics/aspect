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

#ifndef __aspect__vofinterface_QCTrap_h
#define __aspect__vofinterface_QCTrap_h

#include <deal.II/base/config.h>
#include <deal.II/base/quadrature.h>

DEAL_II_NAMESPACE_OPEN

template <int dim>
class QCTrap : public Quadrature<dim>
{
  public:
    QCTrap (const unsigned int n);
};

template <>
QCTrap<1>::QCTrap (const unsigned int n);

DEAL_II_NAMESPACE_CLOSE

#endif
