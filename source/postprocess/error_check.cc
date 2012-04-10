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


#include <aspect/postprocess/error_check.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vectors.h>


#include <math.h>
using namespace dealii;

extern void _Velic_solCx(
    double pos[],
    double _eta_A, double _eta_B,
    double _x_c, int _n,
    double vel[], double* presssure,
    double total_stress[], double strain_rate[] );


namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    class FunctionSolCx : public Function<dim>
    {
      public:
        FunctionSolCx () : Function<dim>() {}

        virtual void vector_value (const Point< dim > & 	p,
                                     Vector< double > & 	values) const
        {
          double pos[2]={p(0),p(1)};
          double result[3];//, vel[2], pressure,
          double total_stress[3], strain_rate[3];
          double eta_A=1.0;
          double eta_B=1e6;

//          printf("t_xx, t_xz fucked !! \n");
//          printf("pos_x, pos_y, vel_x, vel_y, pressure, total_stress[0], total_stress[1], total_stress[2], strain_rate[0], strain_rate[1], strain_rate[2]\n");

              _Velic_solCx(
                pos,
                eta_A, eta_B,
                0.5, 1,
                    &values[0], &values[2], total_stress, strain_rate );

        }
    };


    template <int dim>
    ErrorCheck<dim>::ErrorCheck ()

    {}



    template <int dim>
    std::pair<std::string,std::string>
    ErrorCheck<dim>::execute (TableHandler &statistics)
    {
      FunctionSolCx<dim> ref_func;

      const QGauss<dim> quadrature_formula (this->get_dof_handler().get_fe().base_element(0).degree+2);

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim), dim+2);
      ComponentSelectFunction<dim> comp_p(dim, dim+2);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                             this->get_solution(),
                                             ref_func,
                                             cellwise_errors_u,
                                             quadrature_formula,
                                             VectorTools::L2_norm, &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                             this->get_solution(),
                                             ref_func,
                                             cellwise_errors_p,
                                             quadrature_formula,
                                             VectorTools::L2_norm, &comp_p);

      std::cout << "u_l2:" << cellwise_errors_u.l2_norm() << std::endl;
      std::cout << "p_l2:" << cellwise_errors_p.l2_norm() << std::endl;





              AssertThrow(Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()) == 1,ExcNotImplemented());

              std::ofstream f ((this->get_output_directory() + "output.csv").c_str());
              f.precision (16);
              f << std::scientific;

            //const QGauss<dim> quadrature_formula (this->get_dof_handler().get_fe().base_element(0).degree+2);

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
                            << ' ' << fe_values.JxW (q)
                            << ' ' << velocity_values[q][0]
                            << ' ' << velocity_values[q][1]
                            << ' ' << pressure_values[q]
                            << ' ' << temperature_values[q]
                            << std::endl;
                  }
              }
            return std::make_pair (std::string ("Error Check"), "");
    }


    template <int dim>
    void
    ErrorCheck<dim>::declare_parameters (ParameterHandler &prm)
    {
//      prm.enter_subsection("Postprocess");
//      {
//        prm.enter_subsection("Depth average");
//        {
//          prm.declare_entry ("Time between graphical output", "1e8",
//                             Patterns::Double (0),
//                             "The time interval between each generation of "
//                             "graphical output files. A value of zero indicates "
//                             "that output should be generated in each time step. "
//                             "Units: years if the "
//                             "'Use years in output instead of seconds' parameter is set; "
//                             "seconds otherwise.");
//        }
//        prm.leave_subsection();
//      }
//      prm.leave_subsection();
    }


    template <int dim>
    void
    ErrorCheck<dim>::parse_parameters (ParameterHandler &prm)
    {
//      prm.enter_subsection("Postprocess");
//      {
//        prm.enter_subsection("Depth average");
//        {
//          output_interval = prm.get_double ("Time between graphical output");
//        }
//        prm.leave_subsection();
//      }
//      prm.leave_subsection();
    }





    template <int dim>
    void
    ErrorCheck<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
//      std::ostringstream os;
//      aspect::oarchive oa (os);
//      oa << (*this);

//      status_strings["DepthAverage"] = os.str();
    }


    template <int dim>
    void
    ErrorCheck<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("ErrorCheck") != status_strings.end())
        {

        }

    }



  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ErrorCheck,
                                  "error check",
                                  "A postprocessor that compares the solution to a given analytic solution")
  }
}
