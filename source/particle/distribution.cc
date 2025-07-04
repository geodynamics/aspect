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


namespace aspect
{
  namespace Postprocess
  {

    template <int dim>
    ParticlePDF<dim>::ParticlePDF(const unsigned int granularity,const double bandwidth,KernelFunctions kernel_function)
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
    ParticlePDF<dim>::ParticlePDF(const double bandwidth,KernelFunctions kernel_function)
      :
      bandwidth (bandwidth),
      kernel_function (kernel_function),
      granularity (1),
      max (std::numeric_limits<double>::min()),
      min (std::numeric_limits<double>::max()),
      standard_deviation (numbers::signaling_nan<double>()),
      mean (numbers::signaling_nan<double>()),
      is_defined_per_particle (true)
    {

    }



    template <int dim>
    void ParticlePDF<dim>::fill_from_particle_range(const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range,
                                                    const unsigned int n_particles_in_cell)
    {
      if (is_defined_per_particle == false)
        {
          for (unsigned int x=0; x<this->granularity; ++x)
            {
              for (unsigned int y=0; y<this->granularity; ++y)
                {
                  double granularity_double = static_cast<double>(this->granularity);
                  double x_double = static_cast<double>(x);
                  double y_double = static_cast<double>(y);

                  if (dim == 3)
                    {
                      for (unsigned int z=0; z<this->granularity; ++z)
                        {
                          const double z_double = static_cast<double>(z);
                          Point<dim> reference_point;
                          reference_point[0] = x_double/granularity_double;
                          reference_point[1] = y_double/granularity_double;
                          reference_point[2] = z_double/granularity_double;
                          std::array<unsigned int,dim> table_index;
                          table_index[0] = x;
                          table_index[1] = y;
                          table_index[2] = z;
                          insert_kernel_sum_from_particle_range(reference_point,table_index,n_particles_in_cell,particle_range);
                        }
                    }
                  else
                    {
                      Point<dim> reference_point;
                      reference_point[0] = x_double/granularity_double;
                      reference_point[1] = y_double/granularity_double;
                      std::array<unsigned int,dim> table_index;
                      table_index[0] = x;
                      table_index[1] = y;
                      insert_kernel_sum_from_particle_range(reference_point,table_index,n_particles_in_cell,particle_range);
                    }
                }
            }
        }
      else
        {
          // Sum the value of the kernel function on that position from every other particle.
          for (const auto &reference_particle: particle_range)
            {
              const auto reference_coordinates = reference_particle.get_reference_location();
              double function_value = 0;

              for (const auto &kernel_position_particle: particle_range)
                {
                  const auto kernel_coordinates = kernel_position_particle.get_reference_location();
                  const double distance = reference_coordinates.distance(kernel_coordinates);
                  function_value += apply_selected_kernel_function(distance);
                }

              add_value_to_function_table(function_value/n_particles_in_cell,reference_particle.get_id());
            }

        }


    }

    template <int dim>
    void ParticlePDF<dim>::insert_kernel_sum_from_particle_range(
      const Point<dim> reference_point,
      std::array<unsigned int, dim> table_index,
      const unsigned int n_particles_in_cell,
      const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range)
    {
      for (const auto &particle: particle_range)
        {
          const auto coordinates = particle.get_reference_location();
          const double distance = coordinates.distance(reference_point);
          const double kernel_function_value = apply_selected_kernel_function(distance);
          add_value_to_function_table(table_index,kernel_function_value/n_particles_in_cell);
        }

    }

    template <int dim>
    void ParticlePDF<dim>::add_value_to_function_table(const unsigned int x_index,
                                                       const unsigned int y_index,
                                                       const unsigned int z_index,
                                                       const double input_value)
    {
      Assert(is_defined_per_particle == false, ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));
      TableIndices<dim> entry_index;
      entry_index[0] = x_index;
      entry_index[1] = y_index;
      if (dim == 3)
        entry_index[2] = z_index;

      function_output_table(entry_index) += input_value;
    }



    template <int dim>
    void ParticlePDF<dim>::add_value_to_function_table(const double input_value,types::particle_index reference_particle_id)
    {
      // This needs to be summed before calling this function, unlike
      // when using regular intervals to sum the kernel values.
      Assert(is_defined_per_particle == true, ExcMessage("This function can only be called if the ParticlePDF is computed per particle location."));
      if (input_value > max)
        {
          max = input_value;
          min_particle_index = reference_particle_id;
        }

      if (input_value < min)
        {
          min = input_value;
          min_particle_index = reference_particle_id;
        }
      function_output_vector.push_back(input_value);
    }



    template <int dim>
    void ParticlePDF<dim>::add_value_to_function_table(
      std::array<unsigned int,dim> &index_point,
      const double input_value)
    {
      Assert(is_defined_per_particle == false, ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));
      TableIndices<dim> entry_index;
      entry_index[0] = index_point[0];
      entry_index[1] = index_point[1];
      if (dim == 3)
        entry_index[2] = index_point[2];

      function_output_table(entry_index) += input_value;
    }



    template <int dim>
    double ParticlePDF<dim>::evaluate_function_at_index(
      const unsigned int x_index,
      const unsigned int y_index,
      const unsigned int z_index) const
    {
      Assert(is_defined_per_particle == false, ExcMessage("This function can only be called if the ParticlePDF is computed on a regular grid."));
      TableIndices<dim> entry_index;
      entry_index[0] = x_index;
      entry_index[1] = y_index;
      if (dim == 3)
        entry_index[2] = z_index;

      return function_output_table(entry_index);
    }



    template <int dim>
    void ParticlePDF<dim>::compute_statistical_values()
    {
      standard_deviation = 0;
      mean = 0;

      if (!is_defined_per_particle)
        {
          // Loop through all values of the function to set initial stats
          for (unsigned int x = 0; x< granularity; ++x)
            {
              for (unsigned int y = 0; y< granularity; ++y)
                {
                  if (dim == 3)
                    {
                      for (unsigned int z = 0; z< granularity; ++z)
                        {
                          const double this_value = evaluate_function_at_index(x,y,z);
                          max = std::max(max, this_value);
                          min = std::min(min, this_value);
                          // Sum in mean, then divide after this loop
                          mean += this_value;
                        }
                    }
                  else
                    {
                      const double this_value = evaluate_function_at_index(x,y,0);
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
          for (unsigned int x = 0; x< granularity; ++x)
            {
              for (unsigned int y = 0; y< granularity; ++y)
                {
                  if (dim==3)
                    {
                      for (unsigned int z = 0; z< granularity; ++z)
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
          squared_deviation_sum /= (granularity*dim);
          standard_deviation = std::sqrt(squared_deviation_sum);
        }
      else
        {
          // Loop through the vector and sum mean.
          for (const double this_value : function_output_vector)
            {
              mean += this_value;
            }
          // Set the true mean
          mean /= function_output_vector.size();
          double squared_deviation_sum = 0;
          standard_deviation = 0;

          // Sum all the squared deviations for standard deviation
          for (const double this_value : function_output_vector)
            {
              const double deviation_squared = (this_value-mean)*(this_value-mean);
              squared_deviation_sum += deviation_squared;
            }
          squared_deviation_sum /= function_output_vector.size();
          standard_deviation = std::sqrt(squared_deviation_sum);
        }
    }



    template <int dim>
    double ParticlePDF<dim>::get_max()
    {
      return max;
    }



    template <int dim>
    double ParticlePDF<dim>::get_min()
    {
      return min;
    }



    template <int dim>
    double ParticlePDF<dim>::get_standard_deviation()
    {
      return standard_deviation;
    }



    template <int dim>
    double ParticlePDF<dim>::apply_selected_kernel_function(const double distance) const
    {
      if (kernel_function == KernelFunctions::uniform)
        {
          return kernelfunction_uniform(distance);
        }
      else if (kernel_function == KernelFunctions::triangular)
        {
          return kernelfunction_triangular(distance);
        }
      else if (kernel_function == KernelFunctions::gaussian)
        {
          return kernelfunction_gaussian(distance);
        }
      else if (kernel_function == KernelFunctions::cutoff_function_w1_dealii)
        {
          Functions::CutOffFunctionW1<1> cutoff_function(bandwidth);
          return cutoff_function.value(Point<1>(distance));
        }
      else
        {
          Assert(false, ExcMessage("Unknown kernel function used in insert_kernel_sum_into_pdf."));
          return 0.0;
        }
    }



    template <int dim>
    double ParticlePDF<dim>::kernelfunction_uniform(double distance) const
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
    double ParticlePDF<dim>::kernelfunction_triangular(double distance) const
    {
      if (distance < bandwidth)
        {
          return (1.0-distance)*bandwidth;
        }
      else
        {
          return 0.0;
        }
    }



    template <int dim>
    double ParticlePDF<dim>::kernelfunction_gaussian(double distance) const
    {
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
  namespace Postprocess
  {
#define INSTANTIATE(dim) \
  template class ParticlePDF<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
