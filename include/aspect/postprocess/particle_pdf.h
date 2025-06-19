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

#ifndef _aspect_postprocess_particle_pdf_h
#define _aspect_postprocess_particle_pdf_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <algorithm>
#include <limits>
#include <deal.II/base/table.h>
#include <deal.II/particles/property_pool.h>
namespace aspect
{
  namespace Postprocess
  {
    /**
     * This class handles the point-density function for describing the clustering
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
         * This is the constructor for ParticlePDF.
         * This constructor is called when creating a point-density function which is defined
         * at regular intervals throughout a cell, as opposed to a function which is defined
         * at the positions of each particle in the cell.
         * @param granularity determines the number of reference points to sum the kernel function
         * from. The point density function is only defined at these reference points.
         */
        ParticlePDF(const unsigned int granularity);

        /**
         * This is the constructor for ParticlePDF.
         * This constructor is called when creating a point-density function
         * which is defined at the positions of the particles within the cell,
         * as opposed to at regular intervals based on a granularity value
         */
        ParticlePDF();

        /**
         * Inserts a value into the point-density function.
         * @param x_index The x index in the point-density function to insert at.
         * @param y_index The y index in the point-density function to insert at.
         * @param z_index The z index in the point-density function to insert at.
         * @param input_value The value to insert into the point-density function.
         */
        void add_value_to_function_table(const unsigned int x_index,
                                         const unsigned int y_index,
                                         const unsigned int z_index,
                                         const double input_value);

        /**
         * Inserts a value into the point-density function.
         * @param index_point The index in the point-density function to insert at.
         * @param input_value The value to insert into the point-density function.
         */
        void add_value_to_function_table(
          std::array<unsigned int,dim> &index_point,
          const double input_value);

        /**
         * Inserts a value into the point-density function. This version of the
         * function is intended to be used when `is_defined_per_particle` is `true`.
         * @param input_value The value to insert into the point-density function.
         * @param reference_particle_id the id of the particle from which the input value is taken in reference to.
         */
        void add_value_to_function_table(
          const double input_value,
          types::particle_index reference_particle_id);

        /**
         * Returns the value of the function at the given index.
         * @param x_index The x index in the point-density function to evaluate at.
         * @param y_index The y index in the point-density function to evaluate at.
         * @param z_index The z index in the point-density function to evaluate at.
         */
        double evaluate_function_at_index(
          const unsigned int x_index,
          const unsigned int y_index,
          const unsigned int z_index) const;

        /**
        * Calculates the relevant statistical values from the contents of the PDF.
        */
        void compute_statistical_values();

        /**
         * Returns the maximum of the point-density function.
         */
        double get_max();


        /**
         * Returns the minimum of the point-density function.
         */
        double get_min();

        /**
         * Returns the standard deviation of the point-density function.
         */
        double get_standard_deviation();

      private:
        /**
         * `function_output_table` holds the output of the point-density function
         * at every point where the function is defined, as long as the function
         * has been generated at regular intervals throughout the cell (as opposed
         * to being defined at the position of each particle.
         */
        Table<dim,double> function_output_table;

        /**
         * `pdf_granularity` determines the number of inputs at which
         * the point-density function is defined. The function is defined for
         * pdf_granularity^dim inputs.
         */
        unsigned int pdf_granularity;

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
        std::vector<double> function_output_vector;

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
    };
  }
}

#endif
