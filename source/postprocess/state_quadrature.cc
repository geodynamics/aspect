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


#include <aspect/postprocess/state_quadrature.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <aspect/global.h>
#include <aspect/volume_of_fluid/handler.h>
#include <aspect/volume_of_fluid/utilities.h>
#include <math.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

#include <aspect/simulator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    void
    StateQuadrature<dim>::initialize()
    {
      AssertThrow(dim==2,
                  ExcMessage("Richardson Extrapolation is currently only functional for dim=2."));
      if (refinement_comparison || !amr_comparison)
        {
          AssertThrow((this->get_parameters().initial_adaptive_refinement == 0 ||
                       this->get_parameters().adaptive_refinement_interval == 0),
                      ExcMessage("To compare an AMR solution to the uniform coarser solution "));
        }

      uniform_file_name = this->get_output_directory()
                          + uniform_file_name + "_"
                          + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                          + ".dat";
      coarse_file_name = this->get_output_directory()
                         + coarse_file_name + "_"
                         + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                         + ".dat";
      refined_file_name = this->get_output_directory()
                          + refined_file_name + "_"
                          + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                          + ".dat";
    }

    template <int dim>
    double
    StateQuadrature<dim>::compute_cell_diameter()
    {

      bool initialized = false;
      double min_local_diameter;
      for (const auto &cell: this->get_dof_handler().active_cell_iterators())
        {
          const double diameter = cell->diameter();
          if (!initialized)
            {
              min_local_diameter = diameter;
              initialized = true;
            }
          else
            min_local_diameter = min_local_diameter>diameter?diameter:min_local_diameter;
        }

      const double global_min_diameter = Utilities::MPI::min(min_local_diameter, this->get_mpi_communicator());

      return global_min_diameter;
    }

    template <int dim>
    void
    StateQuadrature<dim>::write_out_header(std::ofstream &ostream)
    {

      // Write header
      ostream << "DIMENSION: " << dim << std::endl;
      ostream << "TIME: " << this->get_time() << std::endl;
      ostream << "CELL DIAMETER: " << cell_diameter << std::endl;
      ostream << "QUADRATURE: " << "Gauss-Legendre:" << quadrature_degree << std::endl;

      // VARIABLE NAMES

      // Write Point header
      for (unsigned int i=0; i< dim; ++i)
        {
          ostream << "X_" << dim_vars[i] << "\t";
        }

      // Write weight header
      // (skip tab because fenceposts at the end of the point values)
      ostream << "W";

      // Write velocity header
      for (unsigned int i=0; i< dim; ++i)
        {
          ostream << "\t" << "V_" << dim_vars[i];
        }

      // Write pressure and temperature to header line
      ostream << "\t" << "P" << "\t" << "T";

      // Write composition vars to header line
      for (auto name: this->introspection().get_composition_names())
        {
          ostream << "\t" << "C_" << name;
        }

      // Write VOF volume fraction to header line
      if (this->get_parameters().volume_of_fluid_tracking_enabled)
        {
          for (unsigned int idx = 0; idx < this->get_volume_of_fluid_handler().get_n_fields(); ++idx)
            {
                ostream << "\t" << "VOF_" << this->get_volume_of_fluid_handler().name_for_field_index(idx);
            }
        }

      ostream << std::endl;
    }

    template <int dim>
    void
    StateQuadrature<dim>::write_out_coarse_data()
    {
      std::ofstream interpolated_data_stream;
      interpolated_data_stream.open(coarse_file_name, std::ios_base::out);
      interpolated_data_stream << std::scientific << std::setprecision(15);

      /**
      * Compute the Legendre gauss points at level 2 indirection.
      **/
      //QGaussLobatto<1> base_quadrature(this->get_stokes_velocity_degree() + 1);
      QGauss<1> base_quadrature(quadrature_degree);
      QIterated<dim> quadrature_rule (base_quadrature, 1);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_rule,
                              update_values |
                              update_quadrature_points |
                              update_JxW_values);

      const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
      const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
      const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;
      
      const unsigned int n_vof_fields = this->get_parameters().volume_of_fluid_tracking_enabled
          ?
          this->get_volume_of_fluid_handler().get_n_fields()
          :
          0;

      write_out_header(interpolated_data_stream);

      // Declaring an iterator over all active cells on local mpi process
      typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active();
      for (; cell != this->get_dof_handler().end();
           ++cell)
        {
          if (!(cell->is_locally_owned()))
            continue;

          fe_values.reinit(cell);
          // Each vector is of the same length.
          const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();

          std::vector<double> interpolated_temperature(quadrature_points.size());
          std::vector<double> interpolated_pressure(quadrature_points.size());
          std::vector<Tensor<1,dim>> interpolated_velocity(quadrature_points.size());
          std::vector<std::vector<double>> interpolated_compositional_fields(this->n_compositional_fields(),
                                                                             std::vector<double>(quadrature_points.size()));
          std::vector<std::vector<double>> interpolated_vof_fields(n_vof_fields,
                                                                   std::vector<double>(quadrature_points.size()));

          fe_values[extractor_pressure].get_function_values(this->get_solution(), interpolated_pressure);
          fe_values[extractor_temperature].get_function_values(this->get_solution(), interpolated_temperature);
          fe_values[extractor_velocity].get_function_values(this->get_solution(), interpolated_velocity);

          typename std::vector<std::vector<double>>::iterator itr_compositional_fields = interpolated_compositional_fields.begin();

          unsigned int index = 0;
          for (; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++)
            {
              fe_values[this->introspection().extractors.compositional_fields[index]].get_function_values(
                this->get_solution(),
                *itr_compositional_fields);
              index++;
            }

          if (this->get_parameters().volume_of_fluid_tracking_enabled)
            {
              for (unsigned int idx = 0; idx < this->get_volume_of_fluid_handler().get_n_fields(); ++idx)
                {
                  const VolumeOfFluidField<dim> &field = this->get_volume_of_fluid_handler().field_struct_for_field_index(idx);
                  const unsigned int volume_of_fluid_component = field.volume_fraction.first_component_index;
                  const FEValuesExtractors::Scalar volume_of_fluid = FEValuesExtractors::Scalar(volume_of_fluid_component);

                  std::vector<double> fraction_values(quadrature_points.size());

                  fe_values[volume_of_fluid].get_function_values(this->get_solution(), fraction_values);

                  const double fraction = fraction_values[0];

                  for (unsigned int i=0; i<quadrature_points.size(); ++i)
                    {
                      interpolated_vof_fields[idx][i] = fraction;
                    }
                }
            }

          typename std::vector<double>::const_iterator itr_temperature = interpolated_temperature.begin();
          typename std::vector<double>::const_iterator itr_pressure = interpolated_pressure.begin();
          typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = interpolated_velocity.begin();

          typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points.begin();

          unsigned int quadrature_point_index = 0;
          for (; itr_quadrature_points != quadrature_points.end();
               itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, quadrature_point_index++)
            {
              // Point
              for (unsigned int i=0; i<dim; ++i)
                {
                  interpolated_data_stream << (*itr_quadrature_points)[i] << "\t";
                }

              // Weight
              interpolated_data_stream << fe_values.JxW(quadrature_point_index);

              for (unsigned int i=0; i<dim; ++i)
                {
                  interpolated_data_stream << "\t" << (*itr_velocity)[i];
                }

              // Pressure and temperature
              interpolated_data_stream << "\t" << *itr_pressure << "\t" << *itr_temperature;

              // Composition fields
              if (this->n_compositional_fields() != 0)
                {
                  unsigned int count = 0;
                  itr_compositional_fields = interpolated_compositional_fields.begin();
                  for (; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++, count++)
                    interpolated_data_stream << "\t" << (*itr_compositional_fields)[quadrature_point_index];
                }

              if (this->get_parameters().volume_of_fluid_tracking_enabled)
                {
                  for (unsigned int i=0; i<this->get_volume_of_fluid_handler().get_n_fields(); ++i)
                    {
                      interpolated_data_stream << "\t" << interpolated_vof_fields[i][quadrature_point_index];
                    }
                }

              interpolated_data_stream << std::endl;
            }
        }

      interpolated_data_stream.close();
    }

    template <int dim>
    void
    StateQuadrature<dim>::write_out_uniform_data()
    {
      std::ofstream interpolated_data_stream;
      interpolated_data_stream.open(uniform_file_name, std::ios_base::out);
      interpolated_data_stream << std::scientific << std::setprecision(15);

      unsigned int uniform_level = this->get_parameters().initial_global_refinement +
                                   this->get_parameters().initial_adaptive_refinement;

      /**
      * Compute the Legendre gauss points at level 2 indirection.
      **/
      //QGaussLobatto<1> base_quadrature(this->get_stokes_velocity_degree() + 1);
      QGauss<1> base_quadrature(quadrature_degree);

      const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
      const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
      const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;
      
      const unsigned int n_vof_fields = this->get_parameters().volume_of_fluid_tracking_enabled
          ?
          this->get_volume_of_fluid_handler().get_n_fields()
          :
          0;

      write_out_header(interpolated_data_stream);

      // Declaring an iterator over all active cells on local mpi process
      typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active();
      for (; cell != this->get_dof_handler().end();
           ++cell)
        {

          if (!(cell->is_locally_owned()))
            continue;

          unsigned int level_diff =  uniform_level - cell->level();
          unsigned int n_divisions = 1;
          for (unsigned int i = 0; i<level_diff; ++i)
            n_divisions *= 2;

          QIterated<dim> quadrature_rule (base_quadrature, n_divisions);

          FEValues<dim> fe_values(this->get_mapping(),
                                  this->get_fe(),
                                  quadrature_rule,
                                  update_values |
                                  update_quadrature_points |
                                  update_JxW_values);

          fe_values.reinit(cell);
          // Each vector is of the same length.
          const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();

          std::vector<double> interpolated_temperature(quadrature_points.size());
          std::vector<double> interpolated_pressure(quadrature_points.size());
          std::vector<Tensor<1,dim>> interpolated_velocity(quadrature_points.size());
          std::vector<std::vector<double>> interpolated_compositional_fields(this->n_compositional_fields(),
                                                                             std::vector<double>(quadrature_points.size()));
          std::vector<std::vector<double>> interpolated_vof_fields(n_vof_fields,
                                                                   std::vector<double>(quadrature_points.size()));

          fe_values[extractor_pressure].get_function_values(this->get_solution(), interpolated_pressure);
          fe_values[extractor_temperature].get_function_values(this->get_solution(), interpolated_temperature);
          fe_values[extractor_velocity].get_function_values(this->get_solution(), interpolated_velocity);

          typename std::vector<std::vector<double>>::iterator itr_compositional_fields = interpolated_compositional_fields.begin();

          unsigned int index = 0;
          for (; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++)
            {
              fe_values[this->introspection().extractors.compositional_fields[index]].get_function_values(
                this->get_solution(),
                *itr_compositional_fields);
              index++;
            }

          if (this->get_parameters().volume_of_fluid_tracking_enabled)
            {
              for (unsigned int idx = 0; idx < this->get_volume_of_fluid_handler().get_n_fields(); ++idx)
                {
                  const VolumeOfFluidField<dim> &field = this->get_volume_of_fluid_handler().field_struct_for_field_index(idx);
                  const unsigned int volume_of_fluid_component = field.volume_fraction.first_component_index;
                  const FEValuesExtractors::Scalar volume_of_fluid = FEValuesExtractors::Scalar(volume_of_fluid_component);

                  std::vector<double> fraction_values(quadrature_points.size());

                  fe_values[volume_of_fluid].get_function_values(this->get_solution(), fraction_values);

                  // We are either at the highest refinement, or not on an interface so use
                  // the volume fraction directly
                  const double fraction = fraction_values[0];

                  for (unsigned int i=0; i<quadrature_points.size(); ++i)
                    {
                      interpolated_vof_fields[idx][i] = fraction;
                    }
                }
            }

          typename std::vector<double>::const_iterator itr_temperature = interpolated_temperature.begin();
          typename std::vector<double>::const_iterator itr_pressure = interpolated_pressure.begin();
          typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = interpolated_velocity.begin();

          typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points.begin();

          unsigned int quadrature_point_index = 0;
          for (; itr_quadrature_points != quadrature_points.end();
               itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, quadrature_point_index++)
            {
              // Point
              for (unsigned int i=0; i<dim; ++i)
                {
                  interpolated_data_stream << (*itr_quadrature_points)[i] << "\t";
                }

              // Weight
              interpolated_data_stream << fe_values.JxW(quadrature_point_index);

              for (unsigned int i=0; i<dim; ++i)
                {
                  interpolated_data_stream << "\t" << (*itr_velocity)[i];
                }

              // Pressure and temperature
              interpolated_data_stream << "\t" << *itr_pressure << "\t" << *itr_temperature;

              // Composition fields
              if (this->n_compositional_fields() != 0)
                {
                  unsigned int count = 0;
                  itr_compositional_fields = interpolated_compositional_fields.begin();
                  for (; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++, count++)
                    interpolated_data_stream << "\t" << (*itr_compositional_fields)[quadrature_point_index];
                }

              if (this->get_parameters().volume_of_fluid_tracking_enabled)
                {
                  for (unsigned int i=0; i<this->get_volume_of_fluid_handler().get_n_fields(); ++i)
                    {
                      interpolated_data_stream << "\t" << interpolated_vof_fields[i][quadrature_point_index];
                    }
                }

              interpolated_data_stream << std::endl;
            }
        }

      interpolated_data_stream.close();
    }

    template <int dim>
    void
    StateQuadrature<dim>::write_out_refined_data()
    {
      std::ofstream interpolated_data_stream;
      interpolated_data_stream.open(refined_file_name, std::ios_base::out);
      interpolated_data_stream << std::scientific << std::setprecision(15);

      /**
      * Compute the Legendre gauss points at level 2 indirection.
      **/
      //QGaussLobatto<1> base_quadrature(this->get_stokes_velocity_degree() + 1);
      QGauss<1> base_quadrature(quadrature_degree);
      QIterated<dim> quadrature_rule (base_quadrature, 2);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_rule,
                              update_values |
                              update_quadrature_points |
                              update_JxW_values);

      const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
      const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
      const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;
      
      std::vector<Point<dim>> vof_refined_unit_recenters(4);

      for (unsigned int i=0; i<2; ++i)
        {
          for (unsigned int j=0; j<2; ++j)
            {
              vof_refined_unit_recenters[2*i+j][0] = -0.25+0.5*i;
              vof_refined_unit_recenters[2*i+j][1] = -0.25+0.5*j;
            }
        }

      std::vector<unsigned int> refined_cell_location(quadrature_rule.size());

      for (unsigned int i=0; i< quadrature_rule.size(); ++i)
        {
          const Point<dim> &point = quadrature_rule.point(i);
          refined_cell_location[i] = 2*((point[0]>0.5)?1:0)+
                                     ((point[1]>0.5)?1:0);
        }
      
      const unsigned int n_vof_fields = this->get_parameters().volume_of_fluid_tracking_enabled
          ?
          this->get_volume_of_fluid_handler().get_n_fields()
          :
          0;

      write_out_header(interpolated_data_stream);

      // Declaring an iterator over all active cells on local mpi process
      typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active();
      for (; cell != this->get_dof_handler().end();
           ++cell)
        {
          if (!(cell->is_locally_owned()))
            continue;

          fe_values.reinit(cell);
          // Each vector is of the same length.
          const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();

          std::vector<double> interpolated_temperature(quadrature_points.size());
          std::vector<double> interpolated_pressure(quadrature_points.size());
          std::vector<Tensor<1,dim>> interpolated_velocity(quadrature_points.size());
          std::vector<std::vector<double>> interpolated_compositional_fields(this->n_compositional_fields(),
                                                                             std::vector<double>(quadrature_points.size()));
          std::vector<std::vector<double>> interpolated_vof_fields(n_vof_fields,
                                                                   std::vector<double>(quadrature_points.size()));

          fe_values[extractor_pressure].get_function_values(this->get_solution(), interpolated_pressure);
          fe_values[extractor_temperature].get_function_values(this->get_solution(), interpolated_temperature);
          fe_values[extractor_velocity].get_function_values(this->get_solution(), interpolated_velocity);

          typename std::vector<std::vector<double>>::iterator itr_compositional_fields = interpolated_compositional_fields.begin();

          unsigned int index = 0;
          for (; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++)
            {
              fe_values[this->introspection().extractors.compositional_fields[index]].get_function_values(
                this->get_solution(),
                *itr_compositional_fields);
              index++;
            }

          if (this->get_parameters().volume_of_fluid_tracking_enabled)
            {
              const std::vector<double> weights = fe_values.get_JxW_values();
              double cell_vol = 0.0;
              for (unsigned int i=0; i<quadrature_points.size(); ++i)
                {
                  cell_vol += weights[i];
                }
              for (unsigned int idx = 0; idx < this->get_volume_of_fluid_handler().get_n_fields(); ++idx)
                {
                  const VolumeOfFluidField<dim> &field = this->get_volume_of_fluid_handler().field_struct_for_field_index(idx);
                  const unsigned int volume_of_fluid_N_component =  field.reconstruction.first_component_index;
                  const FEValuesExtractors::Vector volume_of_fluid_N_n = FEValuesExtractors::Vector(volume_of_fluid_N_component);
                  const FEValuesExtractors::Scalar volume_of_fluid_N_d = FEValuesExtractors::Scalar(volume_of_fluid_N_component+dim);

                  std::vector<Tensor<1,dim> > normal_values(quadrature_points.size());
                  std::vector<double> d_values(quadrature_points.size());

                  fe_values[volume_of_fluid_N_n].get_function_values(this->get_solution(), normal_values);
                  fe_values[volume_of_fluid_N_d].get_function_values(this->get_solution(), d_values);

                  const Tensor<1, dim, double> normal = normal_values[0];
                  const double d = d_values[0];

                  std::vector<double> volume_fractions(vof_refined_unit_recenters.size());

                  for (unsigned int i=0; i<volume_fractions.size(); ++i)
                    {
                      const double d_r = d-normal*vof_refined_unit_recenters[i];

                      volume_fractions[i] = VolumeOfFluid::Utilities::compute_fluid_fraction<dim>(normal, 2.0*d_r);
                    }
                  for (unsigned int i=0; i<quadrature_points.size(); ++i)
                    {
                      interpolated_vof_fields[idx][i] = volume_fractions[refined_cell_location[i]];
                    }
                }
            }

          typename std::vector<double>::const_iterator itr_temperature = interpolated_temperature.begin();
          typename std::vector<double>::const_iterator itr_pressure = interpolated_pressure.begin();
          typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = interpolated_velocity.begin();

          typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points.begin();

          unsigned int quadrature_point_index = 0;
          for (; itr_quadrature_points != quadrature_points.end();
               itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, quadrature_point_index++)
            {
              // Point
              for (unsigned int i=0; i<dim; ++i)
                {
                  interpolated_data_stream << (*itr_quadrature_points)[i] << "\t";
                }

              // Weight
              interpolated_data_stream << fe_values.JxW(quadrature_point_index);

              for (unsigned int i=0; i<dim; ++i)
                {
                  interpolated_data_stream << "\t" << (*itr_velocity)[i];
                }

              // Pressure and temperature
              interpolated_data_stream << "\t" << *itr_pressure << "\t" << *itr_temperature;

              // Composition fields
              if (this->n_compositional_fields() != 0)
                {
                  unsigned int count = 0;
                  itr_compositional_fields = interpolated_compositional_fields.begin();
                  for (; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++, count++)
                    interpolated_data_stream << "\t" << (*itr_compositional_fields)[quadrature_point_index];
                }

              if (this->get_parameters().volume_of_fluid_tracking_enabled)
                {
                  for (unsigned int i=0; i<this->get_volume_of_fluid_handler().get_n_fields(); ++i)
                    {
                      interpolated_data_stream << "\t" << interpolated_vof_fields[i][quadrature_point_index];
                    }
                }

              interpolated_data_stream << std::endl;
            }
        }

      interpolated_data_stream.close();
    }

    template <int dim>
    std::pair<std::string,std::string>
    StateQuadrature<dim>::execute (TableHandler &)
    {
      if ( this->get_time() == end_time)
        {
          cell_diameter = compute_cell_diameter();
          quadrature_degree = this->get_stokes_velocity_degree() + 1;

          if (amr_comparison)
            {
              if ((this->get_parameters().initial_adaptive_refinement == 0 ||
                   this->get_parameters().adaptive_refinement_interval == 0))
                {
                  write_out_coarse_data();

                  return std::pair<std::string, std::string> ("State data written to : ",
                                                              coarse_file_name);
                }
              else
                {
                  write_out_uniform_data();

                  return std::pair<std::string, std::string> ("State data written to : ",
                                                              uniform_file_name);
                }
            }
          else
            {
              if (refinement_comparison)
                {
                  // Write out the current solution, both current and interpolated at a higher resolved mesh.
                  write_out_coarse_data();
                  write_out_refined_data();

                  return std::pair<std::string, std::string> ("State data written to : ",
                                                              coarse_file_name + ", " + refined_file_name);
                }
              else
                {
                  write_out_coarse_data();

                  return std::pair<std::string, std::string> ("State data written to : ",
                                                              coarse_file_name);
                }
            }
        }
      else
        {
          return std::pair<std::string, std::string> ();
        }
    }

    template <int dim>
    void
    StateQuadrature<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("State quadrature");
        {
          prm.declare_entry("File name prefix for interpolated data", "interpolated_data",
                            Patterns::Anything (),
                            "The file prefix for the files containing the"
                            " interpolated solution at the nodal values.");
          prm.declare_entry("AMR comparison", "false",
                            Patterns::Bool (),
                            "Write nodal values to compare to a uniform mesh of"
                            " the same maximum refinement level.");
          prm.declare_entry("Refinement comparison", "false",
                            Patterns::Bool (),
                            "Write nodal values to compare to a uniform mesh of"
                            " one additional refinement level.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    StateQuadrature<dim>::parse_parameters (ParameterHandler &prm)
    {
      bool use_years = prm.get_bool ("Use years in output instead of seconds");

      // Use the existing end time parameter
      //
      // read end time from parameter file. if it is to be interpreted
      // in years rather than seconds, then do the conversion
      end_time = prm.get_double ("End time");
      if ( use_years == true)
        end_time *= year_in_seconds;

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("State quadrature");
        {
          amr_comparison = prm.get_bool("AMR comparison");
          refinement_comparison = prm.get_bool("Refinement comparison");
          std::string file_prefix = prm.get("File name prefix for interpolated data");
          unsigned int max_refinement = this->get_parameters().initial_global_refinement +
                                        this->get_parameters().initial_adaptive_refinement;
          uniform_file_name = file_prefix + "_" + Utilities::int_to_string(max_refinement) + "_a";
          coarse_file_name = file_prefix + "_" + Utilities::int_to_string(max_refinement) + "_c";
          refined_file_name = file_prefix + "_" + Utilities::int_to_string(max_refinement) + "_r";
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(StateQuadrature,
                                  "state quadrature",
                                  "A postprocesser which exports the values of "
                                  "the state variables at the Gauss-Legendre "
                                  "quadrature points. This is useful for "
                                  "postprocessing such as comparing the "
                                  "results to other model configurations as is "
                                  "commonly used for error estimation.")
  }
}
