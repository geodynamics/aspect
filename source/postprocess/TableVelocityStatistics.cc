//-------------------------------------------------------------
//    $Id: velocity_statistics.cc 571 2012-01-03 08:22:42Z geenen $
//
//    Copyright (C) 2011, 2012 by the authors of the ASPECT code
//
//-------------------------------------------------------------

#include <aspect/postprocess/TableVelocityStatistics.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/material_model/table.h>
#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/gravity_model/radial.h>

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
    TableVelocityStatistics<dim>::execute (TableHandler &statistics)
    {
      Assert (dynamic_cast<const MaterialModel::Table<dim> *>(&this->get_material_model())
              != 0,
              ExcMessage ("This postprocessor can only be used when using "
                          "the MaterialModel::Table implementation of the "
                          "material model interface."));
      const MaterialModel::Table<dim> &
      material_model
        = dynamic_cast<const MaterialModel::Table<dim> &>(this->get_material_model());

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

      const double vrms = std::sqrt(global_velocity_square_integral) /
                          std::sqrt(this->get_volume());
      const GeometryModel::SphericalShell<dim> *geometry = dynamic_cast<const GeometryModel::SphericalShell<dim> *>(&this->get_geometry_model());
      const double kappa = material_model.reference_thermal_diffusivity();
      double h;
      if (geometry)
        h = geometry->outer_radius() - geometry->inner_radius();
      else
        h = 1.0; //TODO: define something like geometrymodel::depth?
      double Scaling = h/kappa;
      const double vrmsDimless = Scaling*std::sqrt(global_velocity_square_integral) /
                                 std::sqrt(this->get_volume());
      Scaling = kappa/std::pow(h,2);
      const double timeDimless = Scaling*this->get_time();

      if (this->convert_output_to_years() == true)
        {
          statistics.add_value ("RMS velocity (m/year)",
                                vrms * year_in_seconds);
          statistics.add_value ("Max. velocity (m/year)",
                                global_max_velocity * year_in_seconds);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "RMS velocity (m/year)",
                                      "Max. velocity (m/year)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }
      else
        {
          statistics.add_value ("RMS velocity (m/s)", vrms);
          statistics.add_value ("Max. velocity (m/s)", global_max_velocity);

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "RMS velocity (m/s)",
                                      "Max. velocity (m/s)"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }

      statistics.add_value ("Dimensionless time",
                            timeDimless);
      statistics.set_precision ("Dimensionless time", 6);
      statistics.set_scientific ("Dimensionless time", true);

      statistics.add_value ("Dimensionless RMS velocity",
                            vrmsDimless);
      statistics.set_precision ("Dimensionless RMS velocity", 6);
      statistics.set_scientific ("Dimensionless RMS velocity", true);

      std::ostringstream output;
      output.precision(3);
      if (this->convert_output_to_years() == true)
        output << vrms *year_in_seconds
               << " m/year, "
               << global_max_velocity *year_in_seconds
               << " m/year";
      else
        output << vrms
               << " m/s, "
               << global_max_velocity
               << " m/s";
      if (this->get_time() == 0e0)
        {
          double dT = this->get_boundary_temperature().maximal_temperature() - this->get_boundary_temperature().minimal_temperature();

          Point<dim> representative_point = Point<dim>::unit_vector(dim-1);
          const double gravity = this->get_gravity_model().gravity_vector(representative_point).norm();
          const double Ra = material_model.reference_density()*
                            gravity*
                            material_model.reference_thermal_alpha()*
                            dT*std::pow(h,3)/
                            (material_model.reference_thermal_diffusivity()*
                             material_model.reference_viscosity());

          this->get_pcout()<<  std::endl;
          this->get_pcout()<< "     Reference density (kg/m^3):                    "
                           << material_model.reference_density()
                           << std::endl;
          this->get_pcout()<< "     Reference gravity (m/s^2):                     "
                           << gravity
                           << std::endl;
          this->get_pcout()<< "     Reference thermal expansion (1/K):             "
                           << material_model.reference_thermal_alpha()
                           << std::endl;
          this->get_pcout()<< "     Temperature contrast accross model domain (K): "
                           << dT
                           << std::endl;
          this->get_pcout()<< "     Model domain depth (m):                        "
                           << h
                           << std::endl;
          this->get_pcout()<< "     Reference thermal diffusivity (m^2/s):         "
                           << material_model.reference_thermal_diffusivity()
                           << std::endl;
          this->get_pcout()<< "     Reference viscosity (Pas):                     "
                           << material_model.reference_viscosity()
                           << std::endl;
          this->get_pcout()<< "     Ra number:                                     "
                           << Ra
                           << std::endl;
          this->get_pcout()<<  std::endl;
        }
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
    template class TableVelocityStatistics<deal_II_dimension>;

    ASPECT_REGISTER_POSTPROCESSOR(TableVelocityStatistics,
                                  "velocity statistics for the table model",
                                  "A postprocessor that computes some statistics about the "
                                  "velocity field.")
  }
}
