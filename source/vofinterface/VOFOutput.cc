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

#include <aspect/vofinterface/VOFOutput.h>

namespace aspect
{

  namespace InterfaceTracker
  {
    namespace Output
    {
      using namespace dealii;

      template <int dim>
      VOFOutput<dim>::VOFOutput (const std::string &out_dir,
                                 const MPI_Comm comm)
        : output_dir (out_dir),
          communicator (comm),
          file_index (0),
          times_and_names ()
      {
      }

      template <int dim>
      std::string VOFOutput<dim>::get_filename ()
      {
        const std::string out_file_prefix = ("vofs-"
                                             + Utilities::int_to_string (this->file_index, 5));
        return out_file_prefix + ".vtu";
      }

      template <int dim>
      void VOFOutput<dim>::write_state (const DoFHandler<dim> *dof_handler,
                                        const std::vector<std::string> names,
                                        const std::vector<
                                        DataComponentInterpretation::DataComponentInterpretation> interp,
                                        const Vector<double> &soln,
                                        double time)
      {
        AssertThrow(dof_handler!=NULL, ExcMessage ("DoFHandler must exist."));

        const std::string filename = get_filename ();
        const std::string filepath = this->output_dir + filename;
        const std::string pvd_path = this->output_dir + "vofs.pvd";
        DataOut<dim> data_out;
        DataOutBase::VtkFlags flags;

        flags.time = time;
        flags.cycle = this->file_index;

        data_out.attach_dof_handler (*dof_handler);
        data_out.add_data_vector (soln, names, DataOut<dim>::type_dof_data,
                                  interp);

        data_out.set_flags (flags);
        data_out.build_patches ();

        std::ofstream output (filepath);

        data_out.write_vtu (output);
        ++(this->file_index);

        std::ofstream pvd_output (pvd_path);
        times_and_names.push_back (
          std::pair<double, std::string> (time, filename));
        data_out.write_pvd_record (pvd_output, times_and_names);
      }

      template class VOFOutput<2> ;
      template class VOFOutput<3> ;
    }
  }
}
