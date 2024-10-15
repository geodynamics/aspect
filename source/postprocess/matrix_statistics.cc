/*
  Copyright (C) 2016 - 2023 by the authors of the ASPECT code.

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


#include <aspect/postprocess/matrix_statistics.h>

#include <aspect/simulator.h>
#include <aspect/utilities.h>


namespace
{
  const std::string
  get_stats(const aspect::LinearAlgebra::BlockSparseMatrix &matrix,
            const std::string &matrix_name,
            const MPI_Comm comm)
  {
    std::ostringstream output;

    // convert from bytes into Mb
    const double mb = 1024*1024;
    // sum up local matrix memory usage
    double global_matrix_memory_consumption = dealii::Utilities::MPI::sum(matrix.memory_consumption(),
                                                                          comm);
    output << "\nTotal " << matrix_name << " memory consumption: "
           << std::fixed << std::setprecision(2) << global_matrix_memory_consumption/mb
           << " MB." << std::endl;

    // output number of nonzero elements in matrix. Do so with 1000s separator
    // since they are frequently large; this was previously done by using the empty
    // string locale, but creating std::locale with an empty string caused problems
    // on some platforms, so the functionality to catch the exception and ignore
    // is kept here, even though explicitly setting a facet should always work.
    try
      {
        // Imbue the stream with a locale that does the right thing. The
        // locale is responsible for later deleting the object pointed
        // to by the last argument (the "facet"), see
        // https://en.cppreference.com/w/cpp/locale/locale/locale
        output.imbue(std::locale(std::locale(),
                                 new aspect::Utilities::ThousandSep));
      }
    catch (const std::runtime_error &e)
      {
        // If the locale doesn't work, just give up
      }


    output << "Total " << matrix_name << " nnz: "
           << matrix.n_nonzero_elements() << std::endl;

    // output number of nonzero elements in each matrix block
    output << matrix_name << " nnz by block: " << std::endl;
    for (unsigned int i=0; i<matrix.n_block_rows(); ++i)
      {
        for (unsigned int j=0; j<matrix.n_block_rows(); ++j)
          output << std::setw(12) << matrix.block(i,j).n_nonzero_elements();
        output << std::endl;
      }

    return output.str();
  }
}

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    std::pair<std::string,std::string>
    MatrixStatistics<dim>::execute (TableHandler &)
    {
      std::ostringstream output;
      output << get_stats(this->get_system_matrix(),
                          "system matrix",
                          this->get_mpi_communicator());
      output << get_stats(this->get_system_preconditioner_matrix(),
                          "system preconditioner matrix",
                          this->get_mpi_communicator());

      return std::pair<std::string, std::string> ("",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MatrixStatistics,
                                  "matrix statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the matrices. "
                                  "In particular, it outputs total memory consumption, "
                                  "total non-zero elements, and non-zero elements per "
                                  "block, for system matrix and system preconditioner matrix.")
  }
}
