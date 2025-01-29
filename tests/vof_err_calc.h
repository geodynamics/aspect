/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/interface.h>
#include <aspect/volume_of_fluid/utilities.h>
#include <aspect/volume_of_fluid/handler.h>

// Deal II includes
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    class VolumeOfFluidSpecifiedSolutionDiff : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        VolumeOfFluidSpecifiedSolutionDiff();

        double get_next_t (const double curr_time, const double interval);

        /**
         * Evaluate the solution for some velocity statistics.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        // Error calculation functions
        static std::vector<std::string> error_names ();
        static std::vector<std::string> error_abrev ();
        std::vector<double> calc_error_level_set (const Function<dim> &func,
                                                  const unsigned int n_samp,
                                                  const unsigned int f_ind);

        // Required variables
        bool initialized;
        double error_eval_interval;
        double next_evaluation_time;
        unsigned int n_error_samples;

        /**
         * Level set function for true solution
         */
        std::unique_ptr<Functions::ParsedFunction<dim>> trueSolLS;
    };

    template <int dim>
    VolumeOfFluidSpecifiedSolutionDiff<dim>::VolumeOfFluidSpecifiedSolutionDiff ()
      : initialized(false),
        error_eval_interval (std::numeric_limits<double>::quiet_NaN ()),
        next_evaluation_time (std::numeric_limits<double>::quiet_NaN ()),
        n_error_samples (1)
    {
    }

    template <int dim>
    double VolumeOfFluidSpecifiedSolutionDiff<dim>::get_next_t (const double curr_time,
                                                                const double interval)
    {
      const int i = (int) (curr_time / interval);
      return (i + 1) * interval;
    }

    template <int dim>
    std::pair<std::string, std::string>
    VolumeOfFluidSpecifiedSolutionDiff<dim>::execute (TableHandler &)
    {
      std::string result_string = "";
      std::string label_string = "";
      if (!initialized)
        {
          initialized = true;
          next_evaluation_time = this->get_time ();
        }

      if (this->get_time () >= next_evaluation_time)
        {
          std::vector<std::string> err_abrev = error_abrev ();
          double curr_time = this->get_time();
          if (this->convert_output_to_years())
            curr_time /= year_in_seconds;
          trueSolLS->set_time(curr_time);

          std::ostringstream label_stream;
          std::ostringstream err_str;
          label_stream << "VoF Calculation(";

          unsigned int n_err = err_abrev.size();
          std::vector<std::string>::iterator it = err_abrev.begin();
          for (; it!=err_abrev.end(); )
            {
              label_stream << *it;
              if (++it!=err_abrev.end())
                label_stream << '/';
            }
          label_stream << "):";
          label_string = label_stream.str();

          unsigned int n_volume_of_fluid_fields = this->get_volume_of_fluid_handler().get_n_fields();
          std::vector<double> local_err_vals (n_volume_of_fluid_fields*n_err);
          for (unsigned int f=0; f<n_volume_of_fluid_fields; ++f)
            {
              std::vector<double> l_err_vals_f = calc_error_level_set (*trueSolLS, n_error_samples, f);
              for (unsigned int i=0; i<n_err; ++i)
                local_err_vals[f*n_err+i] = l_err_vals_f[i];
            }

          std::vector<double> global_err_vals(n_volume_of_fluid_fields*n_err);
          Utilities::MPI::sum (local_err_vals, this->get_mpi_communicator(), global_err_vals);

          for (unsigned int f=0; f<n_volume_of_fluid_fields; ++f)
            {
              for (unsigned int i=0; i<n_err; ++i)
                {
                  err_str << std::scientific << std::setprecision (8)
                          << global_err_vals[f*n_err+i];
                  if (i+1<n_err)
                    err_str << " / ";
                }
              if (f+1<n_volume_of_fluid_fields)
                err_str << " // ";
            }
          result_string = err_str.str ();
          next_evaluation_time = get_next_t (this->get_time (), error_eval_interval);
        }
      return std::make_pair (label_string, result_string);
    }

    template <int dim>
    void
    VolumeOfFluidSpecifiedSolutionDiff<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("VoF MMS");
        {
          prm.declare_entry ("Time between error estimates", "1e8",
                             Patterns::Double (0.),
                             "Time interval between error estimates.");

          prm.declare_entry ("Number error samples", "3",
                             Patterns::Integer (1),
                             "Number of iterations for the midpoint quadrature used to accumulate the error estimate");

          prm.enter_subsection ("True LS");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection ();

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    VolumeOfFluidSpecifiedSolutionDiff<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("VoF MMS");
        {
          n_error_samples = prm.get_integer ("Number error samples");

          error_eval_interval = prm.get_double ("Time between error estimates");
          if (this->convert_output_to_years())
            error_eval_interval *= year_in_seconds;
          prm.enter_subsection ("True LS");
          {
            trueSolLS.reset(new Functions::ParsedFunction<dim>(this->get_volume_of_fluid_handler().get_n_fields()));
            trueSolLS->parse_parameters (prm);
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    std::vector<std::string> VolumeOfFluidSpecifiedSolutionDiff<dim>::error_abrev ()
    {
      std::vector<std::string> names (1, "IEstL1");
      names.push_back ("FEstL1");
      return names;
    }

    template <int dim>
    std::vector<double> VolumeOfFluidSpecifiedSolutionDiff<dim>::calc_error_level_set (const Function<dim> &func,
        const unsigned int n_samp,
        const unsigned int f_ind)
    {
      const LinearAlgebra::BlockVector &solution = this->get_solution();

      double I_err_est = 0.0;
      double F_err_est = 0.0;
      const double h = 1.0 / n_samp;

      const DoFHandler<dim> &dof_handler = this->get_dof_handler();
      const FiniteElement<dim> &finite_element = this->get_fe();
      const unsigned int dofs_per_cell = finite_element.dofs_per_cell;
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      const QIterated<dim> quadrature (QMidpoint<1>(), n_samp);
      FEValues<dim> fe_err (this->get_mapping(), finite_element, quadrature,
                            update_JxW_values |
                            update_quadrature_points);

      const unsigned int volume_of_fluid_index = this->get_volume_of_fluid_handler().field_struct_for_field_index(f_ind)
                                                 .volume_fraction.first_component_index;
      const unsigned int volume_of_fluidN_c_index = this->get_volume_of_fluid_handler().field_struct_for_field_index(f_ind)
                                                    .reconstruction.first_component_index;

      // Compare computed solution to reinit using provided function
      for (auto cell : dof_handler.active_cell_iterators ())
        {
          if (!cell->is_locally_owned ())
            continue;

          // Obtain data for this cell
          cell->get_dof_indices (local_dof_indices);
          Tensor<1, dim, double> normal;
          for (unsigned int i=0; i<dim; ++i)
            normal[i] = solution(local_dof_indices[finite_element.component_to_system_index(volume_of_fluidN_c_index+i, 0)]);
          const double d_compute = solution(local_dof_indices[finite_element.component_to_system_index(volume_of_fluidN_c_index+dim, 0)]);
          const double cell_volume_of_fluid = solution(local_dof_indices[finite_element.component_to_system_index(volume_of_fluid_index, 0)]);

          // Calculate approximation for volume
          fe_err.reinit (cell);

          const double cell_vol = cell->measure ();
          const double cell_diam = cell->diameter();
          const double d_true = func.value(cell->barycenter(), f_ind);
          double nnorm1 = 0;
          for (unsigned int i = 0; i < dim; ++i)
            {
              nnorm1 += abs (normal[i]);
            }

          if (abs (d_compute) >= 0.5 * nnorm1 &&
              abs (d_true) >= 0.5 * cell_diam &&
              d_compute * d_true >= 0.0)
            {
              continue;
            }

          double val = 0.0;
          double volume_of_fluid_reinit = 0.0;
          for (unsigned int i = 0; i < fe_err.n_quadrature_points; ++i)
            {
              double d_t = 0.0;
              Tensor<1, dim, double> grad_t;
              Point<dim> xU = quadrature.point (i);
              for (unsigned int di = 0; di < dim; ++di)
                {
                  Point<dim> xH = xU;
                  Point<dim> xL = xU;
                  xH[di] += 0.5 * h;
                  xL[di] -= 0.5 * h;
                  const double dH = func.value(cell->intermediate_point(xH), f_ind);
                  const double dL = func.value(cell->intermediate_point(xL), f_ind);
                  grad_t[di] = (dL - dH);
                  d_t += (0.5 / dim) * (dH + dL);
                }
              const double solution_fluid_fraction_at_point = VolumeOfFluid::Utilities::compute_fluid_fraction (grad_t, d_t);
              volume_of_fluid_reinit += solution_fluid_fraction_at_point*(fe_err.JxW (i) / cell_vol);
              for (unsigned int di = 0; di < dim; ++di)
                {
                  xU[di] -= 0.5;
                }
              const double dot = normal * xU;
              const double computed_fluid_fraction_at_point = VolumeOfFluid::Utilities::compute_fluid_fraction (h * normal,
                                                              (d_compute - dot));
              const double diff = abs (solution_fluid_fraction_at_point - computed_fluid_fraction_at_point);
              val += diff * fe_err.JxW (i);
            }
          I_err_est += val;
          F_err_est += abs (cell_volume_of_fluid - volume_of_fluid_reinit) * cell_vol;
        }

      std::vector<double> errors(1, I_err_est);
      errors.push_back (F_err_est);
      return errors;
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VolumeOfFluidSpecifiedSolutionDiff,
                                  "volume of fluid mms",
                                  "A postprocessor that approximates the interface error.")
  }
}
