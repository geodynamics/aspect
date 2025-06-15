
#ifndef _aspect_postprocess_particle_density_statistics_KDE_h
#define _aspect_postprocess_particle_density_statistics_KDE_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/particles/particle_iterator.h>
#include "particle_density_PDF.h"
#include <deal.II/base/patterns.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the distribution
     * of particles, if possible.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class ParticleDensityStatisticsKDE : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
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
       * DOCUMENTATION NEEDED
       */
      double bandwidth;

      /**
       * The `KernelFunctions` enum class is a data structure which
       * contains the kernel functions available for use in the Kernel
       * Density Estimator.
       */
      enum class KernelFunctions {
        GAUSSIAN,
        EUCLIDEAN,
        TRIANGULAR,
        UNIFORM,
      };

      KernelFunctions kernel_function;

      /**
       * Fills the supplied PDF instance with values from the particles in the given cell.
       * takes a KernelFunction parameter (see enum above). Might use a different method later to pass in different functions.
       */
      void fill_PDF_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell, ParticleDensityPDF<dim> &pdf);
      
      /**
       * This function loops through every particle
       * in the given cell and sums the kernel function from each particle on the given position
       * within the cell. This is only meant to be called from fill_PDF_from_cell. 
       */
      void fill_PDF_point_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell, ParticleDensityPDF<dim> &pdf,const double reference_x,const double reference_y,const double reference_z,const unsigned int table_x,const unsigned int table_y,const unsigned int table_z,const unsigned int particles_in_cell);

      /**
       * sampler X,Y,Z denote the position from which to estimate the kernel function.
       * So as the algorithm iterates through the KDE arrays, samplerX/Y/Z are increased by values according 
       * to the granularity and cell size. There may be a point template class which would be a better input,
       * but for now I am using this and calling the function with samplerZ = 0 for 2D cases.
       */
      double kernelfunction_euclidian(double samplerX, double samplerY, double samplerZ, Particles::ParticleAccessor<dim> particle);

      double kernelfunction_uniform(double samplerX, double samplerY, double samplerZ, Particles::ParticleAccessor<dim> particle);


    };
  }
}


#endif
