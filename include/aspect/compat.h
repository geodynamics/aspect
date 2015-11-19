/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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

#ifndef __aspect__compat_h
#define __aspect__compat_h

#include <aspect/global.h>

/*
 * cross_product
 */
#if !DEAL_II_VERSION_GTE(8,4,0)
template <int dim>
dealii::Tensor<1,dim> cross_product_3d(const dealii::Tensor<1,dim> &a, const dealii::Tensor<1,dim> &b)
{
  Assert (dim==3, dealii::ExcInternalError());
  dealii::Tensor<1,dim> result;
  dealii::cross_product(result, a, b);
  return result;
}

template <int dim>
dealii::Tensor<1,dim> cross_product_2d(const dealii::Tensor<1,dim> &a)
{
  Assert (dim==2, dealii::ExcInternalError());
  dealii::Tensor<1,dim> result;
  dealii::cross_product(result, a);
  return result;
}
#endif


/*
 * MPI::min() functions
 */
#if !DEAL_II_VERSION_GTE(8,3,0)
namespace dealii
{
  namespace Utilities
  {
    namespace MPI
    {
      template <typename T>
      inline
      T min (const T &t,
             const MPI_Comm &mpi_communicator)
      {
        return -max(-t, mpi_communicator);
      }

      template <typename T>
      inline
      void min (const std::vector<T> &values,
                const MPI_Comm       &mpi_communicator,
                std::vector<T>       &minima)
      {
        minima.resize(values.size());
        for (unsigned int i=0; i<values.size(); ++i)
          minima[i] = -values[i];
        max(minima, mpi_communicator, minima);
        for (unsigned int i=0; i<values.size(); ++i)
          minima[i] = -minima[i];
      }
    }
  }
}
#endif

/*
 * Unique_ptr functionality that was introduced in deal.II 8.3 and replaces
 * auto_ptr in CXX17. We need this to silence deprecation warnings in new
 * compilers.
 */
#if !DEAL_II_VERSION_GTE(8,3,0)
#ifdef DEAL_II_WITH_CXX11

#  include <memory>
namespace dealii
{
  namespace std_cxx11
  {
    using std::unique_ptr;
  }
}

#else

#include <boost/scoped_ptr.hpp>

namespace dealii
{
  namespace std_cxx11
  {
    /**
     * Implementation of a basic replacement for C++11's std::unique_ptr class.
     *
     * BOOST does not have a replacement for std::unique_ptr (because unique_ptr
     * requires move semantics that aren't available unless you have a C++11
     * compiler -- in which case you also have std::unique_ptr; see for example
     * http://stackoverflow.com/questions/2953530/unique-ptr-boost-equivalent)
     *
     * Consequently, we emulate the class by just wrapping a boost::scoped_ptr
     * in the cheapest possible way -- by just deriving from it and repeating
     * the basic constructors. Everything else is inherited from the scoped_ptr
     * class.
     *
     * There is no overhead to this approach: scoped_ptr cannot be copied or
     * moved. Instances of unique_ptr cannot be copied, and if you do not have a
     * C++11 compiler, then you cannot move anything anyway.
     */
    template <typename T>
    class unique_ptr : public boost::scoped_ptr<T>
    {
      public:
        unique_ptr () {}

        template<class Y>
        explicit unique_ptr (Y *p)
        :
        boost::scoped_ptr<T>(p)
        {}
    };

  }
}

#endif
#endif



#endif
