//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/velocity_statistics.h>
#include <aspect/simulator.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VelocityStatistics<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim> quadrature_formula (this->get_stokes_dof_handler().get_fe()
                                            .base_element(0).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_stokes_dof_handler().get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<Tensor<1,dim> > velocity_values(n_q_points);

      const FEValuesExtractors::Vector velocities (0);

      double local_velocity_square_integral = 0;
      double local_max_velocity = 0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_stokes_dof_handler().begin_active(),
      endc = this->get_stokes_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[velocities].get_function_values (this->get_stokes_solution(),
                                                       velocity_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                local_velocity_square_integral += ((velocity_values[q] * velocity_values[q]) *
                                                   fe_values.JxW(q));
                local_max_velocity = std::max (std::sqrt(velocity_values[q]*velocity_values[q]),
                                               local_max_velocity);
              }
          }

      const double global_velocity_square_integral
        = Utilities::MPI::sum (local_velocity_square_integral, MPI_COMM_WORLD);
      const double global_max_velocity
        = Utilities::MPI::max (local_max_velocity, MPI_COMM_WORLD);

      const double vrms = std::sqrt(global_velocity_square_integral) / std::sqrt(this->get_volume());
      statistics.add_value ("RMS velocity (cm/year)", vrms * year_in_seconds * 100);
      statistics.add_value ("Max. velocity (cm/year)", global_max_velocity * year_in_seconds * 100);

      std::ostringstream output;
      output.precision(3);
      output << vrms *year_in_seconds * 100
             << " cm/year, "
             << global_max_velocity *year_in_seconds * 100
             << " cm/year";

      return std::pair<std::string, std::string> ("RMS, max velocity:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class VelocityStatistics<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR("velocity statistics", VelocityStatistics)
  }
}
