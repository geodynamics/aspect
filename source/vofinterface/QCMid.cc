/* FEM Interface tracking plugin for ASPECT
 *
 *  Created on: Nov 5, 2015
 *      Author: Jonathan Robey
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

