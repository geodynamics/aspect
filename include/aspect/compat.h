/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#ifndef _aspect_compat_h
#define _aspect_compat_h

#include <aspect/global.h>

// C++11 related includes.
#include <array>
#include <functional>
#include <memory>

#if DEAL_II_VERSION_GTE(9,4,0)
#if !DEAL_II_VERSION_GTE(9,5,0)
#include <aspect/compat/mapping_cartesian.h>
#endif
#endif

namespace big_mpi
{

#if DEAL_II_VERSION_GTE(9,4,0)

  using dealii::Utilities::MPI::broadcast;

#else

  inline MPI_Datatype
  mpi_type_id(const char *)
  {
    return MPI_CHAR;
  }

  /**
   * Broadcast the information in @p buffer from @p root to all
   * other ranks.
   */
  template <typename T>
  void
  broadcast(T                 *buffer,
            const size_t       count,
            const unsigned int root,
            const MPI_Comm    &comm)
  {
    Assert(root < dealii::Utilities::MPI::n_mpi_processes(comm),
           dealii::ExcMessage("Invalid root rank specified."));

    // MPI_Bcast's count is a signed int, so send at most 2^31 in each
    // iteration:
    const size_t max_send_count = std::numeric_limits<signed int>::max();

    size_t total_sent_count = 0;
    while (total_sent_count < count)
      {
        const size_t current_count =
          std::min(count - total_sent_count, max_send_count);

        const int ierr = MPI_Bcast(buffer + total_sent_count,
                                   current_count,
                                   mpi_type_id(buffer),
                                   root,
                                   comm);
        AssertThrowMPI(ierr);
        total_sent_count += current_count;
      }

  }
#endif

}

#endif
