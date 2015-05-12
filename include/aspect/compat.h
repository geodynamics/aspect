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
 * MPI::min() functions
 */
#if DEAL_II_VERSION_GTE(8,3,0)
// use deal.II's MPI::min() functions
#else
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



#endif
