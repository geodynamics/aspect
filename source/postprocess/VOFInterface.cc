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

#include <aspect/postprocess/VOFInterface.h>

namespace aspect
{

  namespace Postprocess
  {
    using namespace dealii;

    template <int dim>
    VOFInterface<dim>::VOFInterface ()
      : voleps (std::numeric_limits<double>::quiet_NaN ()),
        out_interval (std::numeric_limits<double>::quiet_NaN ()),
        next_out_t (0.0),
        n_i_samp (1),
        n_e_samp (1),
        mms (false),
        err_interval (std::numeric_limits<double>::quiet_NaN ()),
        next_err_t (0.0),
        initialized (false),
        output (NULL)
    {
    }

    template <int dim>
    VOFInterface<dim>::~VOFInterface ()
    {
    }

    // Define parameter functions

    template <int dim>
    void VOFInterface<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("VOFInterface");
        {

          prm.declare_entry ("Small volume", "1e-6",
                             Patterns::Double (0, 1),
                             "Minimum significant volume. VOFs below this considered to be zero.");

          prm.declare_entry ("Time between data output", "1e8",
                             Patterns::Double (0), "Time interval between data dumps.");

          prm.declare_entry ("Manufactured Solution", "false",
                             Patterns::Bool (),
                             "When set to true, initial condition function is assumed to be the correct solution for all time and used to calculate the error.");

          prm.declare_entry ("Time between error estimates", "1e8",
                             Patterns::Double (0),
                             "Time interval between error estimates.");

          prm.declare_entry ("Number initialization samples", "3",
                             Patterns::Integer (1),
                             "Number of sampled points per dimension when initializing from VOF");

          prm.declare_entry ("Number error samples", "3",
                             Patterns::Integer (1),
                             "Number of sampled points per dimension when estimating error");

          prm.declare_entry ("Error filename", "",
                             Patterns::FileName (Patterns::FileName::output),
                             "File to write error estimates to for MMS runs, leave blank for no separate record");

          prm.enter_subsection ("Initial conditions");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection ();

        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void VOFInterface<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Postprocess");
      {
        prm.enter_subsection ("VOFInterface");
        {
          voleps = prm.get_double ("Small volume");

          out_interval = prm.get_double ("Time between data output");

          n_i_samp = prm.get_integer ("Number initialization samples");

          n_e_samp = prm.get_integer ("Number error samples");

          mms = prm.get_bool ("Manufactured Solution");

          err_interval = prm.get_double ("Time between error estimates");

          prm.enter_subsection ("Initial conditions");
          {
            initFunc.parse_parameters (prm);
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    double VOFInterface<dim>::get_next_t (double curr_time,
                                          double interval)
    {
      int i = curr_time / interval;
      return (i + 1) * interval;
    }

    //Do execution logic branching

    template <int dim>
    std::pair<std::string, std::string> VOFInterface<dim>::execute (TableHandler &statistics)
    {
      std::string result_string = "";
      if (!initialized)
        {
          // Test implementation requirements
          Assert(dim==2, ExcNotImplemented());

          output = new ::aspect::InterfaceTracker::Output::VOFOutput<dim> (
            this->get_output_directory (), this->get_mpi_communicator ());

          engine.set_voleps (voleps);
          engine.set_triangulation (&(this->get_triangulation ()));
          engine.set_mapping (&(this->get_mapping ()));
          engine.set_parent_dofs (&(this->get_dof_handler ()));
          engine.set_comm (this->get_mpi_communicator ());
          engine.init (initFunc, n_i_samp, this->get_time ());
          result_string += "Initialized.";
          initialized = true;
          next_err_t = this->get_time ();
          next_out_t = this->get_time ();
        }

      if (out_interval > 0 && this->get_time () >= next_out_t)
        {
          engine.calc_normals ();
          result_string += " Wrote " + output->get_filename () + ".";
          output->write_state (engine.get_dof_handler (),
                               ::aspect::InterfaceTracker::VOFEngine<dim>::component_names (),
                               ::aspect::InterfaceTracker::VOFEngine<dim>::component_interpretation (),
                               engine.get_state (), this->get_time ());
          next_out_t = get_next_t (this->get_time (), out_interval);
        }

      if (mms && this->get_time () >= next_err_t)
        {
          std::vector<std::string> err_names = engine.error_names ();
          std::vector<std::string> err_abrev = engine.error_abrev ();
          std::vector<double> err_vals = engine.calc_error (initFunc, n_e_samp,
                                                            this->get_time ());
          std::ostringstream err_str;
          int i=0;
          std::vector<std::string>::iterator it1, it2;
          it2 = err_abrev.begin();
          for (it1=err_names.begin(); it1!=err_names.end(); ++it1,++it2,++i)
            {
              double err = err_vals[i];
              statistics.add_value(*it1, err_vals[i]);
              statistics.set_precision(*it1, 8);
              statistics.set_scientific(*it1, true);
              err_str << " " << *it2 << "=";
              err_str << std::scientific << std::setprecision (15)
                      << err;
            }
          err_str << ".";
          result_string += err_str.str ();
          next_err_t = get_next_t (this->get_time (), err_interval);
        }

      engine.do_step (this->get_solution (), this->get_old_solution (),
                      this->get_timestep ());

      // output->write_state (engine.get_dof_handler (), engine.get_state (),
      //     this->get_time()+this->get_timestep());
      // AssertThrow(false, ExcMessage ("Break for testing"));

      if (result_string != "")
        return std::make_pair ("VOF Calculation:", result_string);
      return std::make_pair ("", "");
    }
  }
}

namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VOFInterface, "vof_interface",
                                  "Postprocessor that tracks fluid interfaces using a VOF method."
                                  "Initial conditions determined by composition fields.")
  }
}
