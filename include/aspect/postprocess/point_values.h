// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_point_values_h
#define _aspect_postprocess_point_values_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that evaluates the solution vector at individual
     * points.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class PointValues : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor
         */
        PointValues ();

        /**
         * Evaluate the solution and determine the values at the
         * selected points.
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
         * Set the time output was supposed to be written. In the simplest
         * case, this is the previous last output time plus the interval, but
         * in general we'd like to ensure that it is the largest supposed
         * output time, which is smaller than the current time, to avoid
         * falling behind with last_output_time and having to catch up once
         * the time step becomes larger. This is done after every output.
         */
        void set_last_output_time (const double current_time);

        /**
         * Interval between the generation of output in seconds.
         *
         * This variable is read from the parameter file through a parameter called 'Time between point values output'.
         */
        double output_interval;

        /**
         * A time (in seconds) the last output has been produced.
         */
        double last_output_time;

        /**
         * Vector of Points representing the points where the solution is to be evaluated
         * that can be used by VectorTools.
         */
        std::vector<Point<dim>> evaluation_points_cartesian;

        /**
         * The values of the solution at the evaluation points.
         */
        std::vector<std::pair<double, std::vector<Vector<double>>>> point_values;

        /**
         * Whether or not to interpret the evaluation points in the input file
         * as natural coordinates or not.
         *
         * This variable is read from the parameter file through a parameter called 'Use natural coordinates'.
         */
        bool use_natural_coordinates;
    };
  }
}


#endif
