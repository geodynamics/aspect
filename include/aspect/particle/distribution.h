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

#ifndef _aspect_particle_distribution_h
#define _aspect_particle_distribution_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <algorithm>
#include <limits>
#include <deal.II/base/table.h>
#include <deal.II/particles/property_pool.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/base/function_lib.h>
#include <vector>

namespace aspect
{
  namespace Particle
  {
    /**
     * This class computes the point-density function that describes the clustering
     * of particles. ParticlePDF holds
     * the values of the point-density function at every point at which
     * the function is defined. The function is defined either at granularity^dim
     * regular intervals throughout each cell, or at the position of
     * each particle in the cell. ParticlePDF computes statistical values
     * from the point-density function including the mean and standard deviation.
     */
    template <int dim>
    class ParticlePDF
    {
      public:

        /**
         * The `KernelFunctions` enum class is a data structure which
         * contains the kernel functions available for use in the Kernel
         * Density Estimator.
         */
        enum class KernelFunction
        {
          gaussian,
          triangular,
          uniform,
          cutoff_function_w1_dealii,
          cutoff_function_c1_dealii
        };

        /**
         * This is the constructor for ParticlePDF.
         * This constructor is called when creating a point-density function which is defined
         * at regular intervals throughout a cell.
         *
         * @param granularity determines the number of reference points at which the point density
         * function is computed. The point density function is only defined at these reference points.
         * @param bandwidth determines the bandwidth to be used in the kernel function.
         * @param kernel_function determines which kernel function to use when generating the point-density function.
         */
        ParticlePDF(const unsigned int granularity,
                    const double bandwidth,
                    const KernelFunction kernel_function);

        /**
         * This is the constructor for ParticlePDF.
         * This constructor is called when creating a point-density function
         * which is defined at the positions of the particles within the cell,
         * as opposed to at regular intervals based on a granularity value.
         * @param bandwidth determines the bandwidth to be used in the kernel function.
         * @param kernel_function determines which kernel function to use when generating the point-density function.
         */
        ParticlePDF(const double bandwidth,
                    const KernelFunction kernel_function);

        /**
         * Fills the point-density function with values from the particles in the given cell.
         * @param particle_range The particle_iterator_range to operate on.
         * @param particle_ranges_to_sum_over The ranges of the particles in neighboring cells.
         * @param n_particles_in_cell The number of particles belonging to the particle manager in question within the cell.
         */
        void
        fill_from_particle_range(const typename Particles::ParticleHandler<dim>::particle_iterator_range particle_range,
                                 const std::vector<typename Particles::ParticleHandler<dim>::particle_iterator_range> particle_ranges_to_sum_over,
                                 const unsigned int n_particles_in_cell);

        /**
         * This function is only called from `fill_from_cell`.
         * It iterates through every particle in the cell and
         * sums the value of the kernel function between the
         * reference point and the position of the cell.
         * @param reference_point The point from which to get the value of the kernel function.
         * @param table_index The index in the PDF to insert the data into.
         * @param n_particles_in_cell The number of particles in the cell.
         * @param particle_range The particle_iterator_range to sum particles from.
         */
        void
        insert_kernel_sum_from_particle_range(const Point<dim> reference_point,
                                              std::array<unsigned int,dim> table_index,
                                              const unsigned int n_particles_in_cell,
                                              const typename Particles::ParticleHandler<dim>::particle_iterator_range particle_range);

        /**
         * Inserts a value into the point-density function.
         * @param x_index The x index in the point-density function to insert at.
         * @param y_index The y index in the point-density function to insert at.
         * @param z_index The z index in the point-density function to insert at.
         * @param input_value The value to insert into the point-density function.
         */
        void
        add_value_to_function_table(const unsigned int x_index,
                                    const unsigned int y_index,
                                    const unsigned int z_index,
                                    const double input_value);

        /**
         * Inserts a value into the point-density function.
         *
         * @param index_point The index in the point-density function to insert at.
         * @param input_value The value to insert into the point-density function.
         */
        void
        add_value_to_function_table(const std::array<unsigned int,dim> &index_point,
                                    const double input_value);

        /**
         * Inserts a value into the point-density function. This version of the
         * function is intended to be used when `is_defined_per_particle` is `true`.
         *
         * @param input_value The value to insert into the point-density function.
         * @param reference_particle_id the id of the particle from which the input value is taken in reference to.
         */
        void
        add_value_to_function_table(const double input_value,
                                    const types::particle_index reference_particle_id);

        /**
         * Returns the value of the function at the given index.
         *
         * @param x_index The x index in the point-density function to evaluate at.
         * @param y_index The y index in the point-density function to evaluate at.
         * @param z_index The z index in the point-density function to evaluate at.
         */
        double
        evaluate_function_at_index(const unsigned int x_index,
                                   const unsigned int y_index,
                                   const unsigned int z_index) const;

        /**
        * Calculates the relevant statistics from the contents of the PDF.
        */
        void
        compute_statistical_values();

        /**
         * Returns the maximum of the point-density function.
         */
        double
        get_max() const;

        /**
         * Returns the minimum of the point-density function.
         */
        double
        get_min() const;

        /**
         * Returns the standard deviation of the point-density function.
         */
        double
        get_standard_deviation() const;

        /**
         * Returns the index of the particle whose position has the highest
         * point-density value. This function only works if the particle density
         * function is defined per particle, instead of being defined on a grid.
         */
        types::particle_index
        get_max_particle() const;

        /**
         * Returns the index of the particle whose position has the lowest
         * point-density value. This function only works if the particle density
         * function is defined per particle, instead of being defined on a grid.
         */
        types::particle_index
        get_min_particle() const;

      private:
        /**
         * `function_output_table` holds the output of the point-density function
         * at every point where the function is defined, as long as the function
         * has been generated at regular intervals throughout the cell (as opposed
         * to being defined at the position of each particle.
         */
        Table<dim,double> function_output_table;

        /**
         * The `bandwidth` variable scales the point-density function.
         * Choosing an appropriate bandwidth is important because a
         * bandwidth value which is too low or too high can result
         * in oversmoothing or undersmoothing of the point-density function.
         * Oversmoothing or undersmoothing results in a function which
         * represents the underlying data less accurately.
         */
        double bandwidth;

        /**
         * `kernel_function` is an internal variable to keep track of which
         * kernel function is being used by the ParticlePDF.
         */
        KernelFunction kernel_function;

        /**
         * `granularity` determines the number of inputs at which
         * the point-density function is defined. The function is defined for
         * granularity^dim inputs.
         */
        unsigned int granularity;

        /**
         * `max` holds the maximum value of the point-density function after it has been computed.
         */
        double max;

        /**
         * `min` holds the minimum value of the point-density function after it has been computed.
         */
        double min;

        /**
         * `standard_deviation` holds the standard_deviation of the point-density function after it has been computed.
         */
        double standard_deviation;

        /**
         * `mean` holds the mean value of the point-density function after it has been computed.
         */
        double mean;

        /**
         * `function_output_vector` holds the output generated by summing the kernel functions
         * from the positions of each particle in the cell, as opposed to summing the function
         * at regular intervals throughout the cell.
         */
        small_vector<double> function_output_vector;

        /**
         * `is_defined_per_particle` keeps track of whether or not this point-density function
         * is defined at regular intervals throughout the cell or if the function is defined at
         * the position of each particle in the cell.
         */
        bool is_defined_per_particle;

        /**
         * The index of the particle in the cell with the maximum
         * point-density function value.
         */
        types::particle_index max_particle_index;

        /**
         * The index of the particle in the cell with the minimum
         * point-density function value.
         */
        types::particle_index min_particle_index;

        /**
         * Returns the value of the selected kernel function.
         *
         * @param distance the distance to pass to the selected kernel function.
         */
        double apply_selected_kernel_function(const double distance) const;

        /**
         * The Uniform kernel function returns a value of 1.0 as long as the
         * distance is less than the KDE's bandwidth.
         *
         * @param distance the output of the kernel function depends on the distance between the reference point and the center of the kernel function.
         */
        double kernelfunction_uniform(const double distance) const;

        /**
         * The Triangular kernel function returns a value of 1.0 minus
         * the distance variable.
         *
         * @param distance the output of the kernel function depends on the distance between the reference point and the center of the kernel function.
         */
        double kernelfunction_triangular(const double distance) const;

        /**
         * The gaussian function returns the value of a gaussian distribution
         * at the specified distance.
         *
         * @param distance the output of the kernel function depends on the distance between the reference point and the center of the kernel function.
         */
        double kernelfunction_gaussian(const double distance) const;
    };
  }
}

#endif
