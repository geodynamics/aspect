/*
 * FEM VOF Interface tracking plugin for Aspect
 *
 * Copyright (C) 2015 Jonathan M Robey
 *
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
