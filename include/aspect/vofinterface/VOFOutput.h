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

#ifndef __aspect__vofinterface_VOFOutput_h
#define __aspect__vofinterface_VOFOutput_h

// Aspect includes
#include <aspect/global.h>

// Deal II includes
#include <deal.II/base/data_out_base.h>

#include <deal.II/numerics/data_out.h>

// Definition

namespace aspect
{
  namespace InterfaceTracker
  {
    namespace Output
    {
      using namespace dealii;

      template <int dim>
      class VOFOutput
      {
        private:
          const std::string output_dir;
          const MPI_Comm communicator;
          unsigned int file_index;
          std::vector<std::pair<double, std::string>> times_and_names;

        public:
          VOFOutput (const std::string &out_dir,
                     const MPI_Comm comm);

          std::string get_filename ();

          void write_state (const DoFHandler<dim> *handler,
                            const std::vector<std::string> names,
                            const std::vector<
                            DataComponentInterpretation::DataComponentInterpretation> interp,
                            const Vector<double> &soln,
                            double time);
      };
    }
  }
}

#endif
