/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */

#include <aspect/material_model/tan_gurnis.h>
#include <aspect/postprocess/tan_gurnis.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    std::pair<std::string,std::string>
    TanGurnis<dim>::execute (TableHandler &statistics)
    {
      AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,
                  ExcNotImplemented());

      const MaterialModel::TanGurnis<dim> *
      material_model = dynamic_cast<const MaterialModel::TanGurnis<dim> *>(&this->get_material_model());

      AssertThrow(material_model!=NULL, ExcMessage("tan gurnis postprocessor only works with tan gurnis material model"));

      double ref=1.0/this->get_triangulation().begin_active()->minimum_vertex_distance();
      std::ofstream f ((this->get_output_directory() + "vel_" +
                        Utilities::int_to_string(static_cast<unsigned int>(ref)) +
                        ".csv").c_str());
      f.precision (16);
      f << material_model->parameter_Di() << ' '
        << material_model->parameter_gamma() << ' '
        << material_model->parameter_wavenumber() << ' '
        << material_model->parameter_a() << ' '
        << " -1 -1 -1" << std::endl; //pad to 7 values, so matlab is happy

      f << std::scientific;


      const QGauss<dim> quadrature_formula (this->get_fe().base_element(0).degree+2);

      const unsigned int n_q_points =  quadrature_formula.size();
      FEValues<dim> fe_values (this->get_mapping(), this->get_fe(),  quadrature_formula,
                               update_JxW_values | update_values | update_quadrature_points);

      std::vector<Tensor<1, dim> > velocity_values (quadrature_formula.size());
      std::vector<double>         temperature_values (quadrature_formula.size());
      std::vector<double>         pressure_values (quadrature_formula.size());

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell != endc; ++cell)
        {
          fe_values.reinit (cell);
          fe_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(), velocity_values);
          fe_values[this->introspection().extractors.pressure].get_function_values (this->get_solution(), pressure_values);
          fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(), temperature_values);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              f
                  <<  fe_values.quadrature_point (q) (0)
                  << ' ' << fe_values.quadrature_point (q) (1)
                  << ' ' << velocity_values[q][0]
                  << ' ' << velocity_values[q][1]
                  << ' ' << fe_values.JxW (q)
                  << ' ' << pressure_values[q]
                  << ' ' << temperature_values[q]
                  << std::endl;
            }
        }

      return std::make_pair("writing:", "output.csv");
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(TanGurnis,
                                  "Tan Gurnis error",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Tan/Gurnis (2007) paper with the one computed by ASPECT "
                                  "by outputing data that is compared using a matlab script.")
  }
}
