// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_stokes_residual_h
#define _aspect_postprocess_stokes_residual_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>
#include <deal.II/lac/solver_control.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that outputs Stokes residuals to a file stokes_residuals.txt.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class StokesResidual : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        StokesResidual ();

        /**
         * This attaches to the Stokes solver signal.
         */
        void initialize() override;

        /**
         * Generate the output file.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

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
         * A structure for a single time step record.
         */
        struct DataPoint
        {
          double time;
          unsigned int solve_index;
          std::vector<double> values;

          template <class Archive>
          void serialize (Archive &ar, const unsigned int version);
        };

        /**
         * Callback function to collect the data.
         */
        void stokes_solver_callback (const SolverControl &solver_control_cheap,
                                     const SolverControl &solver_control_expensive);

        /**
         * An array of all the past values
         */
        std::vector<DataPoint> entries;

    };
  }
}


#endif
