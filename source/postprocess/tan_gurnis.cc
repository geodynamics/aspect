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

//#include <aspect/material_model/duretz_et_al.h>
#include <aspect/postprocess/tan_gurnis.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vectors.h>

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

      std::ofstream f ((this->get_output_directory() + "output.csv").c_str());
                    f.precision (16);
                    f << "0 0 0 0"
//                         EquationData::MaterialModel::Di << ' '
//                       << EquationData::MaterialModel::gamma << ' '
//                       << EquationData::MaterialModel::wavenumber << ' '
//                       << EquationData::MaterialModel::a
                       << " -1 -1 -1" << std::endl; //pad to 7 values, so matlab is happy

                    f << std::scientific;


                  const QGauss<dim> quadrature_formula (this->get_dof_handler().get_fe().base_element(0).degree+2);

                  const unsigned int n_q_points =  quadrature_formula.size();
                  FEValues<dim> fe_values (this->get_mapping(), this->get_dof_handler().get_fe(),  quadrature_formula,
                                           update_JxW_values | update_values | update_quadrature_points);

                  const unsigned int dofs_per_cell = fe_values.get_fe().dofs_per_cell;
                  const FEValuesExtractors::Vector velocities (0);
                  const FEValuesExtractors::Scalar pressure (dim);
                  const FEValuesExtractors::Scalar temperature (dim+1);

                  std::vector<Tensor<1, dim> > velocity_values (quadrature_formula.size());
                  std::vector<double>         temperature_values (quadrature_formula.size());
                  std::vector<double>         pressure_values (quadrature_formula.size());

                  typename DoFHandler<dim>::active_cell_iterator
                  cell = this->get_dof_handler().begin_active(),
                  endc = this->get_dof_handler().end();
                  for (; cell != endc; ++cell)
                    {
                      fe_values.reinit (cell);
                      fe_values[velocities].get_function_values (this->get_solution(), velocity_values);
                      fe_values[pressure].get_function_values (this->get_solution(), pressure_values);
                      fe_values[temperature].get_function_values (this->get_solution(), temperature_values);

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
                                  "by outputing data that is compared using a matlab script.");
  }
}

