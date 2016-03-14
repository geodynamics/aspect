/*
 * FEM VOF Interface tracking plugin for Aspect
 *
 * Copyright (C) 2015 Jonathan M Robey
 *
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
