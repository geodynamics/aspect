//-------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/temperature_statistics.h>
#include <aspect/simulator.h>

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
    TemperatureStatistics<dim>::execute (TableHandler &statistics)
    {
//TODO: think about whether it would be useful to only get the degree of the temperature component of the FESystem
      const QGauss<dim> quadrature_formula (this->get_dof_handler().get_fe().degree+1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_dof_handler().get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      const FEValuesExtractors::Scalar temperature (dim+1);
      std::vector<double> temperature_values(n_q_points);

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

      double local_temperature_integral = 0;

      // compute the integral quantities by quadrature
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[temperature].get_function_values (this->get_solution(),
							temperature_values);
            for (unsigned int q=0; q<n_q_points; ++q)
              local_temperature_integral += temperature_values[q]*fe_values.JxW(q);
          }

      // compute min/max by simply
      // looping over the elements of the
      // solution vector. the reason is
      // that minimum and maximum are
      // usually attained at the
      // boundary, and so taking their
      // values at Gauss quadrature
      // points gives an inaccurate
      // picture of their true values
      double local_min_temperature = std::numeric_limits<double>::max();
      double local_max_temperature = std::numeric_limits<double>::min();
      for (unsigned int i=0; i<this->get_solution().block(2).local_size(); ++i)
        {
          local_min_temperature
            = std::min<double> (local_min_temperature,
                                this->get_solution().block(2).trilinos_vector()[0][i]);
          local_max_temperature
            = std::max<double> (local_max_temperature,
                                this->get_solution().block(2).trilinos_vector()[0][i]);
        }

      const double global_temperature_integral
        = Utilities::MPI::sum (local_temperature_integral, MPI_COMM_WORLD);
      double global_min_temperature = 0;
      double global_max_temperature = 0;

      // now do the reductions that are
      // min/max operations. do them in
      // one communication by multiplying
      // one value by -1
      {
        double local_values[2] = { -local_min_temperature, local_max_temperature };
        double global_values[2];

        Utilities::MPI::max (local_values, MPI_COMM_WORLD, global_values);

        global_min_temperature = -global_values[0];
        global_max_temperature = global_values[1];
      }

      statistics.add_value ("Minimal temperature (K)",
                            global_min_temperature);
      statistics.add_value ("Average temperature (K)",
                            global_temperature_integral / this->get_volume());
      statistics.add_value ("Maximal temperature (K)",
                            global_max_temperature);
      statistics.add_value ("Average nondimensional temperature (K)",
                            global_temperature_integral / this->get_volume() /
                            this->get_boundary_temperature().maximal_temperature());

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Minimal temperature (K)",
                                  "Average temperature (K)",
                                  "Maximal temperature (K)",
                                  "Average nondimensional temperature (K)"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << global_min_temperature << " K, "
             << global_temperature_integral / this->get_volume() << " K, "
             << global_max_temperature << " K";

      return std::pair<std::string, std::string> ("Temperature min/avg/max:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    template class TemperatureStatistics<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR(TemperatureStatistics,
                                  "temperature statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the temperature field.")
  }
}
