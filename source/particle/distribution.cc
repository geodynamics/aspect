/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/particle/distribution.h>
#include <aspect/particle/manager.h>
#include <deal.II/fe/mapping.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim>
    ParticlePDF<dim>::ParticlePDF(const unsigned int granularity,
                                  const double bandwidth,
                                  KernelFunction kernel_function)
      :
      bandwidth (bandwidth),
      kernel_function (kernel_function),
      granularity (granularity),
      max (std::numeric_limits<double>::min()),
      min (std::numeric_limits<double>::max()),
      standard_deviation (numbers::signaling_nan<double>()),
      mean (numbers::signaling_nan<double>()),
      is_defined_per_particle (false)
    {
      TableIndices<dim> bucket_sizes;
      for (unsigned int i=0; i<dim; ++i)
        {
          bucket_sizes[i] = granularity;
        }
      function_output_table.reinit(bucket_sizes);
    }



    template <int dim>
    ParticlePDF<dim>::ParticlePDF(const double bandwidth,
                                  KernelFunction kernel_function)
      :
      bandwidth (bandwidth),
      kernel_function (kernel_function),
      granularity (1),
      max (std::numeric_limits<double>::min()),
      min (std::numeric_limits<double>::max()),
      standard_deviation (numbers::signaling_nan<double>()),
      mean (numbers::signaling_nan<double>()),
      is_defined_per_particle (true)
    {}



    template <int dim>
    void
    ParticlePDF<dim>::fill_from_particle_range(const typename Particles::ParticleHandler<dim>::particle_iterator_range &particle_range,
                                               const std::vector<typename Particles::ParticleHandler<dim>::particle_iterator_range> &particle_ranges_to_sum_over,
                                               const unsigned int n_particles_in_cell,
                                               const typename dealii::Mapping<dim> &mapping,
                                               const typename Triangulation<dim>::active_cell_iterator &cell)
    {
      if (is_defined_per_particle == false)
        {
          for (unsigned int x=0; x<granularity; ++x)
            {
              for (unsigned int y=0; y<granularity; ++y)
                {
                  const double spacing = 1./granularity;

                  if (dim == 3)
                    {
                      for (unsigned int z=0; z<granularity; ++z)
                        {
                          const Point<dim> reference_point = Point<dim>(x*spacing + spacing/2., y*spacing + spacing/2., z*spacing + spacing/2.);
                          std::array<unsigned int,dim> table_index;
                          table_index[0] = x;
                          table_index[1] = y;
                          table_index[2] = z;
                          insert_kernel_sum_from_particle_range(reference_point,table_index,cell,particle_ranges_to_sum_over,mapping);
                        }
                    }
                  else
                    {
                      const Point<dim> reference_point = Point<dim>(x*spacing + spacing/2., y*spacing + spacing/2.);
                      std::array<unsigned int,dim> table_index;
                      table_index[0] = x;
                      table_index[1] = y;
                      insert_kernel_sum_from_particle_range(reference_point,table_index,cell,particle_ranges_to_sum_over,mapping);
                    }
                }
            }
        }
      else if (is_defined_per_particle == true)
        {
          // Sum the value of the kernel function on that position from every other particle.
          for (const auto &reference_particle: particle_range)
            {
              // Cell diameter is used to normalize the distance value by the size of the cell containing the
              // particle whose PDF value is being calculated.
              const double cell_diameter_scaled_to_dimensions = cell->diameter() / std::sqrt(dim);
              const auto &particle_coordinates = reference_particle.get_location();
              double function_value = 0;

              for (const auto &particle_range_to_sum: particle_ranges_to_sum_over)
                {
                  for (const auto &kernel_position_particle: particle_range_to_sum)
                    {
                      const auto &kernel_coordinates = kernel_position_particle.get_location();
                      const double distance = particle_coordinates.distance(kernel_coordinates);
                      const double distance_normalized = distance/cell_diameter_scaled_to_dimensions;
                      function_value += apply_selected_kernel_function(distance_normalized);
                    }
                }

              add_value_to_function_table(function_value/n_particles_in_cell,reference_particle.get_id());
            }

        }


    }



    template <int dim>
    void
    ParticlePDF<dim>::insert_kernel_sum_from_particle_range(const Point<dim> &reference_point,
                                                            const std::array<unsigned int, dim> &table_index,
                                                            const typename Triangulation<dim>::active_cell_iterator &cell,
                                                            const std::vector<typename Particles::ParticleHandler<dim>::particle_iterator_range> &particle_ranges_to_sum_over,
                                                            const typename dealii::Mapping<dim> &mapping)
    {
      Assert(is_defined_per_particle == false,
             ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));

      TableIndices<dim> entry_index;
      entry_index[0] = table_index[0];
      entry_index[1] = table_index[1];
      if (dim == 3)
        entry_index[2] = table_index[2];
      function_output_table(entry_index) = 0.0;
      unsigned int n_particles_in_ranges = 0;

      // The reference point is in coordinates relative to the cell, we need it in real coordinates.
      const Point<dim> reference_point_real_coordinates = mapping.transform_unit_to_real_cell(cell,reference_point);

      for (const auto &particle_range_to_sum: particle_ranges_to_sum_over)
        {
          for (const auto &particle: particle_range_to_sum)
            {
              ++n_particles_in_ranges;
              const double cell_diameter_scaled_to_dimensions = cell->diameter() / std::sqrt(dim);
              const double distance = reference_point_real_coordinates.distance(particle.get_location());
              const double distance_normalized = distance/cell_diameter_scaled_to_dimensions;
              const double kernel_function_value = apply_selected_kernel_function(distance_normalized);
              add_value_to_function_table(table_index,kernel_function_value);
            }
        }

      if (n_particles_in_ranges>0)
        {
          // Scale KDE value by number of particle contributions
          function_output_table(entry_index) /= static_cast<double>(n_particles_in_ranges);
        }
      // else: do nothing as the initial PDF value of 0.0 is correct

    }



    template <int dim>
    void
    ParticlePDF<dim>::add_value_to_function_table(const unsigned int x_index,
                                                  const unsigned int y_index,
                                                  const unsigned int z_index,
                                                  const double input_value)
    {
      Assert(is_defined_per_particle == false,
             ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));

      TableIndices<dim> entry_index;
      entry_index[0] = x_index;
      entry_index[1] = y_index;
      if (dim == 3)
        entry_index[2] = z_index;

      function_output_table(entry_index) += input_value;
    }



    template <int dim>
    void
    ParticlePDF<dim>::add_value_to_function_table(const double input_value,
                                                  const types::particle_index reference_particle_id)
    {
      Assert(is_defined_per_particle == true,
             ExcMessage("This function can only be called if the ParticlePDF is computed per particle location."));

      if (input_value > max)
        {
          max = input_value;
          max_particle_index = reference_particle_id;
        }
      if (input_value < min)
        {
          min = input_value;
          min_particle_index = reference_particle_id;
        }

      function_output_vector.push_back(input_value);
    }



    template <int dim>
    void
    ParticlePDF<dim>::add_value_to_function_table(const std::array<unsigned int,dim> &index_point,
                                                  const double input_value)
    {
      Assert(is_defined_per_particle == false,
             ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));

      TableIndices<dim> entry_index;
      entry_index[0] = index_point[0];
      entry_index[1] = index_point[1];
      if (dim == 3)
        entry_index[2] = index_point[2];

      function_output_table(entry_index) += input_value;
    }



    template <int dim>
    double
    ParticlePDF<dim>::evaluate_function_at_index(const unsigned int x_index,
                                                 const unsigned int y_index,
                                                 const unsigned int z_index) const
    {
      Assert(is_defined_per_particle == false,
             ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));

      TableIndices<dim> entry_index;
      entry_index[0] = x_index;
      entry_index[1] = y_index;
      if (dim == 3)
        entry_index[2] = z_index;

      return function_output_table(entry_index);
    }



    template <int dim>
    void
    ParticlePDF<dim>::compute_statistical_values()
    {
      standard_deviation = 0.0;
      mean = 0.0;

      if (!is_defined_per_particle)
        {
          const double spacing = 1./granularity;

          // Loop through all values of the function to set initial stats
          for (unsigned int x=0; x<granularity; ++x)
            {
              for (unsigned int y=0; y<granularity; ++y)
                {
                  if (dim == 3)
                    {
                      for (unsigned int z=0; z<granularity; ++z)
                        {
                          const double this_value = evaluate_function_at_index(x,y,z);

                          const Point<dim> position_in_cell = Point<dim>(x*spacing + spacing/2., y*spacing + spacing/2., z*spacing +spacing/2);

                          // Record the positions of max and min values as well. These are useful for adding particles.
                          if (this_value > max)
                            {
                              max_position = position_in_cell;
                            }

                          // if this value is less than the min, we found a new minimum value, clear the list
                          if (this_value < min)
                            {
                              min_positions.clear();
                              min_positions.push_back(position_in_cell);
                            }
                          else if (this_value == min)
                            {
                              // Add this position to the vector of identical minimum values
                              min_positions.push_back(position_in_cell);
                            }

                          max = std::max(max, this_value);
                          min = std::min(min, this_value);
                          // Sum in mean, then divide after this loop
                          mean += this_value;
                        }
                    }
                  else
                    {
                      const double this_value = evaluate_function_at_index(x,y,0);
                      const Point<dim> position_in_cell = Point<dim>(x*spacing + spacing/2., y*spacing + spacing/2.);

                      //record the positions of max and min values as well. These are useful for adding particles.
                      if (this_value > max)
                        {
                          max_position = position_in_cell;
                        }
                      // if this value is less than the min, we found a new minimum value, clear the list
                      if (this_value < min)
                        {
                          min_positions.clear();
                          min_positions.push_back(position_in_cell);
                        }
                      else if (this_value == min)
                        {
                          // Add this position to the vector of identical minimum values
                          min_positions.push_back(position_in_cell);
                        }

                      max = std::max(max, this_value);
                      min = std::min(min, this_value);

                      // Sum in mean, then divide after this loop
                      mean += this_value;
                    }
                }
            }
          // Set the true mean
          mean /= Utilities::fixed_power<dim>(granularity);
          double squared_deviation_sum = 0;

          // Then sum all the squared deviations for standard deviation
          for (unsigned int x=0; x<granularity; ++x)
            {
              for (unsigned int y=0; y<granularity; ++y)
                {
                  if (dim == 3)
                    {
                      for (unsigned int z=0; z<granularity; ++z)
                        {
                          TableIndices<dim> entry_index;
                          entry_index[0] = x;
                          entry_index[1] = y;
                          entry_index[2] = z;
                          const double this_value = function_output_table(entry_index);
                          const double deviation_squared = (this_value-mean)*(this_value-mean);
                          squared_deviation_sum += deviation_squared;
                        }
                    }
                  else
                    {
                      TableIndices<dim> entry_index;
                      entry_index[0] = x;
                      entry_index[1] = y;
                      const double this_value = function_output_table(entry_index);
                      const double deviation_squared = (this_value-mean)*(this_value-mean);
                      squared_deviation_sum += deviation_squared;
                    }
                }
            }

          // Standard deviation of all the defined points in the density function
          const double squared_deviation_mean = squared_deviation_sum / Utilities::fixed_power<dim>(granularity);
          standard_deviation = std::sqrt(squared_deviation_mean);
        }
      else
        {
          // Loop through the vector and compute the sum.
          for (const double this_value : function_output_vector)
            mean += this_value;

          // Compute the mean
          mean /= function_output_vector.size();
          double squared_deviation_sum = 0;

          // Sum all the squared deviations for standard deviation
          for (const double this_value : function_output_vector)
            {
              const double deviation_squared = (this_value-mean)*(this_value-mean);
              squared_deviation_sum += deviation_squared;
            }
          const double squared_deviation_mean = squared_deviation_sum / function_output_vector.size();
          standard_deviation = std::sqrt(squared_deviation_mean);
        }
    }



    template <int dim>
    double
    ParticlePDF<dim>::get_max() const
    {
      return max;
    }



    template <int dim>
    Point<dim>
    ParticlePDF<dim>::get_max_position() const
    {
      return max_position;
    }



    template <int dim>
    std::vector<Point<dim>>
    ParticlePDF<dim>::get_min_positions() const
    {
      return min_positions;
    }



    template <int dim>
    double
    ParticlePDF<dim>::get_min() const
    {
      return min;
    }



    template <int dim>
    double
    ParticlePDF<dim>::get_standard_deviation() const
    {
      return standard_deviation;
    }



    template <int dim>
    types::particle_index ParticlePDF<dim>::get_max_particle() const
    {
      return max_particle_index;
    }



    template <int dim>
    types::particle_index ParticlePDF<dim>::get_min_particle() const
    {
      return min_particle_index;
    }



    template <int dim>
    double
    ParticlePDF<dim>::apply_selected_kernel_function(const double distance) const
    {
      if (kernel_function == KernelFunction::uniform)
        {
          return kernelfunction_uniform(distance);
        }
      else if (kernel_function == KernelFunction::triangular)
        {
          return kernelfunction_triangular(distance);
        }
      else if (kernel_function == KernelFunction::gaussian)
        {
          return kernelfunction_gaussian(distance);
        }
      else if (kernel_function == KernelFunction::cutoff_function_w1_dealii)
        {
          Functions::CutOffFunctionW1<1> cutoff_function(bandwidth);
          return cutoff_function.value(Point<1>(distance));
        }
      else if (kernel_function == KernelFunction::cutoff_function_c1_dealii)
        {
          Functions::CutOffFunctionC1<1> cutoff_function(bandwidth);
          return cutoff_function.value(Point<1>(distance));
        }
      else
        {
          Assert(false, ExcMessage("Unknown kernel function used in apply_selected_kernel_function."));
          return 0.0;
        }
    }



    template <int dim>
    double
    ParticlePDF<dim>::kernelfunction_uniform(double distance) const
    {
      if (distance < bandwidth)
        {
          return 1;
        }
      else
        {
          return 0.0;
        }
    }



    template <int dim>
    double
    ParticlePDF<dim>::kernelfunction_triangular(double distance) const
    {
      if (distance < bandwidth)
        {
          return (bandwidth-distance)/bandwidth;
        }
      else
        {
          return 0.0;
        }
    }



    template <int dim>
    double
    ParticlePDF<dim>::kernelfunction_gaussian(double distance) const
    {
      // This implementation can be problematic. It ensures the function
      // has a compact support (within bandwidth), but it includes a jump
      // when distance == bandwidth, since a gaussian usually extends to
      // infinity
      if (distance < bandwidth)
        {
          const double exponent = distance * distance / (2.*bandwidth*bandwidth);
          const double gaussian = (1 / (bandwidth * std::sqrt(2*numbers::PI))) * std::exp(-exponent);
          return gaussian;
        }
      else
        {
          return 0.0;
        }
    }
  }
}

// explicit instantiation
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class ParticlePDF<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
