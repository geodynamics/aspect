/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_mobility_statistics_h
#define _aspect_postprocess_mobility_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about mobility following Tackley (2000) and Lourenco et al. (2020)
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class MobilityStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for some mobility statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;
 
        double get_average_mobility() const;
        
        double get_combined_mobility() const;     
     
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
         
        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

       private:
         /**
          * Set the time output was supposed to be written. In the simplest
          * case, this is the previous last output time plus the interval, but
          * in general we'd like to ensure that it is the largest supposed
          * output time, which is smaller than the current time, to avoid
          * falling behind with last_output_time and having to catch up once
          * the time step becomes larger. This is done after every output.
          */
          void set_last_output_time (const double current_time);
          void set_last_average_time (const double current_time);
    
          /**
           * Interval between the generation of output in seconds.
           */
          double output_interval;
          double average_interval;
     
          /**
           * A time (in seconds) the last output has been produced.
           */
          double last_output_time;
          double last_average_time;
 
          double combined_mobility;

          double average_mobility; 
    };
  }
}


#endif
