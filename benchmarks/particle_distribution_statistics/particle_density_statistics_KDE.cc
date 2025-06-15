#include "particle_density_statistics_KDE.h"
#include <aspect/particle/manager.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/tria_accessor.h>
#include <cmath>

namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    ParticleDensityStatisticsKDE<dim>::execute (TableHandler &statistics)
    {
      unsigned int cells_with_particles = 0;
      double standard_deviation_sum = 0;
      double standard_deviation_min = std::numeric_limits<double>::max();
      double standard_deviation_max = std::numeric_limits<double>::min();
      double function_min_sum = 0;
      double function_max_sum = 0;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
      {
        if (cell->is_locally_owned())
        {
          unsigned int particles_in_cell = 0;
          for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
            particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
          
          //call fill_PDF_from_cell 1 time per cell
          if(particles_in_cell > 0)
          {
            cells_with_particles++;
            ParticleDensityPDF pdf = ParticleDensityPDF<dim>(granularity);
            fill_PDF_from_cell(cell,pdf);
            pdf.set_statistical_values();
            if (pdf.standard_deviation > standard_deviation_max)
              standard_deviation_max = pdf.standard_deviation;
            if (pdf.standard_deviation < standard_deviation_min)
              standard_deviation_min = pdf.standard_deviation;

            standard_deviation_sum += pdf.standard_deviation;
            function_min_sum += pdf.min;
            function_max_sum += pdf.max;
          }
        }
      }
     
      // Get final values from all processors
      const double global_standard_deviation_max = Utilities::MPI::max (standard_deviation_max, this->get_mpi_communicator());
      const double global_standard_deviation_min = Utilities::MPI::min (standard_deviation_min, this->get_mpi_communicator());
      const double global_cells_with_particles = Utilities::MPI::sum (cells_with_particles, this->get_mpi_communicator());
      const double global_standard_deviation_sum = Utilities::MPI::sum (standard_deviation_sum, this->get_mpi_communicator());
      const double global_standard_deviation_mean = global_standard_deviation_sum/global_cells_with_particles;
  
      // Get the average of the functions max and min values
      const double global_function_min_sum = Utilities::MPI::sum (function_min_sum, this->get_mpi_communicator());
      const double global_function_max_sum = Utilities::MPI::sum (function_max_sum, this->get_mpi_communicator());
      const double global_function_min_mean = global_function_min_sum/global_cells_with_particles;
      const double global_function_max_mean = global_function_max_sum/global_cells_with_particles;

      // write to statistics file
      statistics.add_value ("Minimum PDF standard deviation ", global_standard_deviation_min);
      statistics.add_value ("Mean of PDF standard deviation: ", global_standard_deviation_mean);
      statistics.add_value ("Maximum PDF standard deviation: ", global_standard_deviation_max);
      statistics.add_value ("Mean of PDF minimum values: ", global_function_min_mean);
      statistics.add_value ("Maximum PDF maximum values: ", global_function_max_mean);


      std::ostringstream output;
      output << global_standard_deviation_min <<"," << global_standard_deviation_mean <<"," << global_standard_deviation_max <<",";
      output << "\n" << "Particle Distribution Statistics (mean of function min/max): " << global_function_min_mean <<"," << global_function_max_mean;


      return std::pair<std::string, std::string> ("Particle Distribution Statistics (function standard deviation min/mean/max):",
                                                  output.str());
    }



    template <int dim>
    std::list<std::string>
    ParticleDensityStatisticsKDE<dim>::required_other_postprocessors() const
    {
      return {"particles"};
    }
    


    template <int dim>
    void ParticleDensityStatisticsKDE<dim>::fill_PDF_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell,ParticleDensityPDF<dim> &pdf)
    {
      unsigned int particles_in_cell = 0;
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
      {
        particles_in_cell += this->get_particle_manager(particle_manager_index).get_particle_handler().n_particles_in_cell(cell);
      }

      for(unsigned int x=0;x<granularity;x++)
      {
        for(unsigned int y=0;y<granularity;y++)
        {
          const double reference_x = x/granularity;
          const double reference_y = y/granularity;
          if (dim == 3) {
            for(unsigned int z=0;z<granularity;z++)
            {
              const double reference_z = z/granularity;
              fill_PDF_point_from_cell(cell,pdf,reference_x,reference_y,reference_z,x,y,z,particles_in_cell);
            }
          } else {
            fill_PDF_point_from_cell(cell,pdf,reference_x,reference_y,0,x,y,0,particles_in_cell);
          }
        }
      }
    }



    template <int dim>
    void ParticleDensityStatisticsKDE<dim>::fill_PDF_point_from_cell(const typename Triangulation<dim>::active_cell_iterator &cell, ParticleDensityPDF<dim> &pdf, const double reference_x, const double reference_y, const double reference_z,const unsigned int table_x,const unsigned int table_y,const unsigned int table_z,const unsigned int particles_in_cell)
    {
      //in here, add every particle in the cell using the kernel function and add that value to the PDF
      for (unsigned int particle_manager_index = 0; particle_manager_index < this->n_particle_managers(); ++particle_manager_index)
      {
        //only use the particles in the cell by using particles_in_cell(cell)
        auto first_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).begin();
        auto last_particle_in_cell = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell).end();

        const typename Particle::ParticleHandler<dim>::particle_iterator_range particle_range = this->get_particle_manager(particle_manager_index).get_particle_handler().particles_in_cell(cell);

        for (const auto &particle: particle_range) {
          if (kernel_function == KernelFunctions::UNIFORM){
            const double PDF_value = kernelfunction_uniform(reference_x,reference_y,reference_z,particle);
            pdf.add_value_to_function_table(table_x,table_y,table_z,PDF_value/particles_in_cell);
          } 
          else if (kernel_function == KernelFunctions::EUCLIDEAN) {//gaussian not implemented yet.
            const double PDF_value = kernelfunction_euclidian(reference_x,reference_y,reference_z,particle);
            pdf.add_value_to_function_table(table_x,table_y,table_z,PDF_value/particles_in_cell);
          }          
          else if (kernel_function == KernelFunctions::GAUSSIAN) {//gaussian not implemented yet.
           // const double PDF_value = kernelfunction_euclidian(reference_x,reference_y,reference_z,particle);
            //pdf.add_value_to_function_table(table_x,table_y,table_z,PDF_value/particles_in_cell);
          } else { //default to euclidean
            //const double PDF_value = kernelfunction_euclidian(reference_x,reference_y,reference_z,particle);
            //pdf.add_value_to_function_table(table_x,table_y,table_z,PDF_value/particles_in_cell);
          }
        }
      }
    }


    // this function is called from getPDF, if getPDF is called with the KernelFunctions::Euclidean parameter.
    // in 2d, this would be a circle? in 3d a sphere?
    template <int dim>
    double ParticleDensityStatisticsKDE<dim>::kernelfunction_euclidian(double samplerX, 
                                                                       double samplerY, 
                                                                       double samplerZ, 
                                                                       Particles::ParticleAccessor<dim> particle)
    {
        const auto coordinates = particle.get_reference_location();
        const double particle_x = coordinates[0];
        const double particle_y = coordinates[1];
        const double particle_z = coordinates[1];
        const double distance = std::sqrt(((samplerX-particle_x)*(samplerX-particle_x))+((samplerY-particle_y)*(samplerY-particle_y))+((samplerZ-particle_z)*(samplerZ-particle_z)));
        const double distance_times_bandwith = distance*bandwidth;
        return distance_times_bandwith;
    }

    template <int dim>
    double ParticleDensityStatisticsKDE<dim>::kernelfunction_uniform(double samplerX, 
                                                                     double samplerY, 
                                                                     double samplerZ, 
                                                                     Particles::ParticleAccessor<dim> particle)
    {
        const auto coordinates = particle.get_reference_location();
        const double particle_x = coordinates[0];
        const double particle_y = coordinates[1];
        const double particle_z = coordinates[1];
        //const double distance = ((samplerX-particle_x)*(samplerX-particle_x))+((samplerY-particle_y)*(samplerY-particle_y))+((samplerZ-particle_z)*(samplerZ-particle_z));
       
        // I think this should be equivalent to a uniform function in 3d.
        // https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
        bool sampler_within_function = true;
        if (std::abs(samplerX-particle_x)>bandwidth){
          sampler_within_function = false;
        }
        if (std::abs(samplerY-particle_y)>bandwidth){
          sampler_within_function = false;
        }    
        if (dim==3 && (std::abs(samplerZ-particle_z))>bandwidth){
          sampler_within_function = false;
        }
      
        // With a uniform function, anything within the function outputs a value of 0.5. According to wikipedia...
        if (sampler_within_function == true){
          return 0.5;
        } else {
          return 0.0;
        } 
    }

    template <int dim>
    void
    ParticleDensityStatisticsKDE<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle Density KDE");
        {
         prm.declare_entry("Kernel Function","Uniform",
                            Patterns::Selection("Uniform|Gaussian|Triangular|Euclidean"),
                            "The kernel smoothing function to use for kernel density estimation."
                            );
          prm.declare_entry("KDE Granularity","2",
                            Patterns::Integer (1),
                            "The granularity parameter determines how many discrete inputs exist for "
                            "the probability density function generated by the kernel density estimator. "
                            "The domain of the function is multidimensional so the granularity value determines "
                            "the range of inputs in each dimension. For example, a granularity value of 2 "
                            "results in a PDF which is defined for the inputs 0-1 in each of its dimensions. ");
           prm.declare_entry("Kernel Bandwidth","0.3",
                            Patterns::Double (0.01),
                            "The granularity parameter determines how many discrete inputs exist for "
                            "the probability density function generated by the kernel density estimator. "
                            "The domain of the function is multidimensional so the granularity value determines "
                            "the range of inputs in each dimension. For example, a granularity value of 2 "
                            "results in a PDF which is defined for the inputs 0-1 in each of its dimensions. ");
     
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ParticleDensityStatisticsKDE<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Particle Density KDE");
        {
          //KDE_per_particle = prm.get_integer("Use KDE Per Particle");
          granularity = prm.get_integer("KDE Granularity");
          bandwidth = prm.get_double("Kernel Bandwidth");
          std::string kernel_function_string = prm.get("Kernel Function");

          if (kernel_function_string =="Triangular")
          {
            kernel_function = KernelFunctions::TRIANGULAR;
          } 
          else if (kernel_function_string =="Gaussian")
          {
            kernel_function = KernelFunctions::GAUSSIAN;
          } 
          else if (kernel_function_string =="Euclidean")
          {
            kernel_function = KernelFunctions::EUCLIDEAN;
          }           
          else
          {
            kernel_function = KernelFunctions::UNIFORM;
          }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(ParticleDensityStatisticsKDE,
                                  "Particle Density KDE",
                                  "A postprocessor that computes some statistics about "
                                  "the particle distribution, if present in this simulation. ")
  }
}
