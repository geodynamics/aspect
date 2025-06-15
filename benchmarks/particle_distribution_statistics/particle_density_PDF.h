#ifndef _aspect_postprocess_particle_density_PDF_h
#define _aspect_postprocess_particle_density_PDF_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/particles/particle_iterator.h>
#include <vector>
#include <array>
#include <map>//I guess later I could try unordered_map which might be faster.
#include <stdexcept>
#include <deal.II/base/table.h>

namespace aspect
{
  namespace Postprocess
  {
    /**
     * A class to handle holding the KDE function. eventually should be templated so it can handle 3D too.
     *
     */
    template <int dim>
    class ParticleDensityPDF
    {
      public:
        Table<dim,double> function_output_table;
        unsigned int pdf_granularity;
        double max,min,standard_deviation,mean,median,quartile_first,quartile_second,quartile_IQR;



        ParticleDensityPDF(unsigned int granularity)
        {
            this->pdf_granularity = granularity;
            max = std::numeric_limits<double>::min();
            min = std::numeric_limits<double>::max();
            standard_deviation = -1; //-1 is not really possible, so if this value is returned, we know something is wrong
            mean = -1;
            median = -1;
            quartile_first = -1;
            quartile_second = -1;
            quartile_IQR = -1;

            TableIndices<dim> bucket_sizes;
            for (unsigned int i=0; i<dim; ++i){
              bucket_sizes[i] = granularity;
            }
            function_output_table.reinit(bucket_sizes);
        };



        void add_value_to_function_table(const unsigned int x_index, const unsigned int y_index, const unsigned int z_index,const double input_value)
        {
          TableIndices<dim> entry_index;
          entry_index[0] = x_index;
          entry_index[1] = y_index;
          if (dim == 3)
            entry_index[2] = z_index;

          function_output_table(entry_index) += input_value;
        };



        double evaluate_function_at_index(const unsigned int x_index, const unsigned int y_index, const unsigned int z_index) const
        {
          TableIndices<dim> entry_index;
          entry_index[0] = x_index;
          entry_index[1] = y_index;
          if (dim == 3)
            entry_index[2] = z_index;

          return function_output_table(entry_index);
        }



        /*
        setStatisticalValues just sets max, min, std deviation, mean, median, q1,q2,iqr all at once. 
        needs to be called once the pdf is filled though.
        */
        void set_statistical_values()
        {
          max = std::numeric_limits<double>::min();
          min = std::numeric_limits<double>::max();
          standard_deviation = 0;
          mean = 0;
          median = 0;
          quartile_first = 0;
          quartile_second = 0;
          quartile_IQR = 0;

          // loop through all values of the function to set initial stats.
          for(unsigned int x = 0; x< pdf_granularity;x++)
          {
            for(unsigned int y = 0; y< pdf_granularity;y++)
            {
              if (dim == 3){
                for(unsigned int z = 0; z< pdf_granularity;z++)
                {
                  double this_value = evaluate_function_at_index(x,y,z);
                  if (this_value > max)
                    max = this_value;
      
                  if (this_value < min)
                    min = this_value;

                  //sum in mean, then divide after this loop
                  mean += this_value;
                }
              } else {
                  double this_value = evaluate_function_at_index(x,y,0);

                  max = std::max(max, this_value);
                  min = std::min(min, this_value);

                  //sum in mean, then divide after this loop
                  mean += this_value;
              }
            }
          }
          // set the true mean
          mean /= (Utilities::fixed_power<dim>(pdf_granularity));//think this should be the total number of points.
          double squared_deviation_sum = 0;

          // this sum all the squared deviations for standard deviation.
          for(unsigned int x = 0; x< pdf_granularity;x++)
          {
            for(unsigned int y = 0; y< pdf_granularity;y++)
            {
              TableIndices<dim> entry_index;
              entry_index[0] = x;
              entry_index[1] = y;
              double this_value = function_output_table(entry_index);
              double deviation_squared = (this_value-mean)*(this_value-mean);
              squared_deviation_sum += deviation_squared;
            }
          }
          // standard deviation of all the defined points in the density function
          squared_deviation_sum /= (pdf_granularity*dim);
          standard_deviation = std::sqrt(squared_deviation_sum);
        };  



   };
  }
}

#endif
