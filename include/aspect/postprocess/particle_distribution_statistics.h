
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

#ifndef _aspect_postprocess_particle_distribution_statistics_h
#define _aspect_postprocess_particle_distribution_statistics_h

#include <aspect/particle/distribution.h>
#include <deal.II/base/function_lib.h>
#include <aspect/particle/manager.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that generates a point-density function describing
     * the distribution of particles within cells. The point-density function
     * is generated using a Kernel Density Estimator. The postprocessor calculates
     * standard deviation and max/min of the point-density functions per cell.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class ParticleDistributionStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize function. Called once at the beginning of the model.
         */
        void initialize() override;

        /**
         * Evaluate the solution for some particle statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Let the postprocessor manager know about the other postprocessors
         * this one depends on. Specifically, the particles postprocessor.
         */
        std::list<std::string>
        required_other_postprocessors() const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * If `KDE_per_particle` is true, the point-density function is defined
         * at the position of every particle in the cell. If it is false, the
         * point-density-function is defined on a regular grid throughout the cell.
         */
        bool KDE_per_particle;

        /**
         * The `granularity` variable determines the number of points at which
         * the point-density function is defined within each cell. For example,
         * a value of 2 means that the point-density function is defined at
         * $2\times 2=4$ points in 2D. This variable only applies if
         * `KDE_per_particle` is false.
         */
        unsigned int granularity;

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
         * kernel function was read from the .prm file.
         */
        typename Particle::ParticlePDF<dim>::KernelFunction kernel_function;

        /**
         * This function returns a vector containing the particle_iterator_ranges of
         * the cells neighboring the supplied cell. This is used to sum the
         * kernel function using particles across cell boundaries.
         * @param cell The cell to find neighboring cell's particles for.
         * @param particle_handler The particle handler of the current particle manager
         */
        std::vector<typename Particles::ParticleHandler<dim>::particle_iterator_range>
        get_neighboring_particle_ranges(const typename Triangulation<dim>::active_cell_iterator &cell,
                                        const typename Particles::ParticleHandler<dim> &particle_handler);

        /**
         * Fills the supplied PDF instance with values from the particles in the given cell.
         *
         * @param cell The cell whose particles to populate the PDF with.
         * @param pdf The instance of ParticlePDF to fill with data from the cell's particles.
         */
        void
        fill_PDF_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell,
                           Particle::ParticlePDF<dim> &pdf);

        /**
         * This function is only called from `fill_PDF_from_cell`.
         * It iterates through every particle in the cell and
         * sums the value of the kernel function between the
         * reference point and the position of the cell.
         *
         * @param cell The cell whose particles to populate the PDF with.
         * @param reference_point The point from which to get the value of the kernel function.
         * @param table_index The index in the PDF to insert the data into.
         * @param particles_in_cell The number of particles in the cell.
         * @param pdf The instance of ParticlePDF to fill with data from the cell's particles.
         */
        void insert_kernel_sum_into_pdf(const typename Triangulation<dim>::active_cell_iterator &cell,
                                        const Point<dim> reference_point,
                                        const std::array<unsigned int,dim> table_index,
                                        const unsigned int particles_in_cell,
                                        Particle::ParticlePDF<dim> &pdf);

        /**
         * Returns the value of the selected kernel function.
         *
         * @param distance the distance to pass to the selected kernel function.
         */
        double
        apply_selected_kernel_function(const double distance) const;

        /**
         * The Uniform kernel function returns a value of 1.0 as long as the
         * distance is less than the KDE's bandwidth.
         *
         * @param distance the output of the kernel function depends on the distance between the reference point and the center of the kernel function.
         */
        double
        kernelfunction_uniform(const double distance) const;

        /**
         * The Triangular kernel function returns a value of 1.0 minus
         * the distance variable.
         *
         * @param distance the output of the kernel function depends on the distance between the reference point and the center of the kernel function.
         */
        double
        kernelfunction_triangular(const double distance) const;

        /**
         * The gaussian function returns the value of a gaussian distribution
         * at the specified distance.
         *
         * @param distance the output of the kernel function depends on the distance between the reference point and the center of the kernel function.
         */
        double
        kernelfunction_gaussian(const double distance) const;

        /**
         * Cached information that stores information about the grid so that we
         * do not need to recompute it every time properties_at_points() is called.
         */
        std::unique_ptr<GridTools::Cache<dim>> grid_cache;
    };
  }
}


#endif
