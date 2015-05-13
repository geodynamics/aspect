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


namespace dealii
{
  namespace Utilities
  {
    namespace MPI
    {
      namespace
      {

        inline MPI_Datatype mpi_type_id (const unsigned int *)
        {
          return MPI_UNSIGNED;
        }


        inline MPI_Datatype mpi_type_id (const unsigned long int *)
        {
          return MPI_UNSIGNED_LONG;
        }
      }

      template <typename T>
      inline
      T min (const T &t,
             const MPI_Comm &mpi_communicator)
      {
        T result;
        MPI_Allreduce (const_cast<void *>(static_cast<const void *>(&t)),
                       &result, 1, mpi_type_id(&t), MPI_MIN,
                       mpi_communicator);
        return result;
      }
    }
  }
}

#endif
