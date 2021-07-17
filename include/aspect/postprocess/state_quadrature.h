/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__postprocess_state_quadrature_h
#define __aspect__postprocess_state_quadrature_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

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
    class StateQuadrature : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        void initialize ();

        /**
         * Evaluate the solution and determine the values at the
         * selected points.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        virtual
        double
        compute_cell_diameter();

        /**
         * Write ascii format data-file header.
         */
        void
        write_out_header(std::ofstream &ostream);

        /**
         * Write out the solution data at current mesh scale. Again, this is
         * written out into an ascii file and instead, should be written out in
         * binary format.
         */
        void
        write_out_coarse_data();

        /**
         * Write out the solution data at current mesh scale for quadrature
         * points corresponding to a uniform mesh. Again, this is written out
         * into an ascii file and instead, should be written out in binary
         * format.
         */
        void
        write_out_uniform_data();

        /**
         * Write out the solution data at mesh scale for one more level of
         * refinement. Again, this is written out into an ascii file and
         * instead, should be written out in binary format.
         */
        void
        write_out_refined_data();

      private:
        /**
         * Run time parameters.
         */
        const std::string dim_vars[3] = {"X", "Y", "Z"};
        std::string uniform_file_name;
        std::string coarse_file_name;
        std::string refined_file_name;
        unsigned int quadrature_degree;
        double cell_diameter;
        double end_time;
        bool amr_comparison;
        bool refinement_comparison;
    };
  }
}


#endif
