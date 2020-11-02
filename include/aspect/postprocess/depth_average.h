/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_depth_average_h
#define _aspect_postprocess_depth_average_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that generates depth average output in periodic
     * intervals or every time step.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class DepthAverage : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        DepthAverage ();

        /**
         * Evaluate the solution and compute the requested depth averages.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

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
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * Interval between the generation of output in seconds. This parameter is read
         * from the input file and consequently is not part of the state that
         * needs to be saved and restored.
         */
        double output_interval;

        /**
         * A time (in seconds) the last output has been produced.
         */
        double last_output_time;

        /**
         * The formats in which to produce graphical output. This also
         * determines the extension of the file names to which to write.
         */
        std::vector<std::string> output_formats;

        /**
         * Number of zones in depth direction over which we are supposed to
         * average.
         */
        unsigned int n_depth_zones;

        /**
         * The boundaries of the depth zones.
         * This vector contains exactly n_depth_zones + 1 entries,
         * and entries i and i+1 represent the lower and upper depth bound
         * of zone i.
         */
        std::vector<double> depth_bounds;

        /**
         * List of the quantities to calculate for each depth zone.
         */
        std::vector<std::string> variables;

        /**
         * A structure for a single time step record.
         */
        struct DataPoint
        {
          double time;
          std::vector<std::vector<double> > values;

          template <class Archive>
          void serialize (Archive &ar, const unsigned int version);
        };

        /**
         * An array of all the past values
         */
        std::vector<DataPoint> entries;

        /**
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);
    };
  }
}


#endif
