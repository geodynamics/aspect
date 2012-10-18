/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: zero_velocity.h 1071 2012-08-01 16:50:02Z bangerth $  */


#ifndef __aspect__velocity_boundary_conditions_gplates_h
#define __aspect__velocity_boundary_conditions_gplates_h

#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator.h>


namespace aspect
{
  namespace VelocityBoundaryConditions
  {
    using namespace dealii;

    namespace internal
    {
      class GPlatesLookup;
    }

    /**
     * A class that implements prescribed velocity boundary conditions
     * determined from a GPlates input files.
     *
     * @ingroup VelocityBoundaryConditionsModels
     */
    template <int dim>
    class GPlates : public Interface<dim>
    {
      public:
        /**
         * Constructor.
         */
        GPlates ();

        /**
         * Return the boundary velocity as a function of position. For the
         * current class, this function returns value from gplates.
         */
        Tensor<1,dim>
        boundary_velocity (const Point<dim> &position) const;

        /**
         * Initialization function. This function is called once at the beginning
         * of the program and loads the first set of GPlates files describing initial
         * conditions at the start time.
         */
        void
        initialize (const GeometryModel::Interface<dim> &geometry_model);

        /**
        * A function that is called at the beginning of each time step
        * to indicate what the model time is for which the boundary
        * values will next be evaluated. Also loads the next velocity
        * files if necessary and outputs a warning if the end of the set of
        * velocity files if reached.
        */
        void
        set_current_time (const double time);

        /**
        * Declare the parameters this class takes through input files.
        */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter
         * file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Pointer to the geometry object in use.
         */
        const GeometryModel::Interface<dim> *geometry_model;

        /**
        * A variable that stores the current time of the simulation. Derived
        * classes can query this variable. It is set at the beginning of each
        * time step.
        */
        double current_time;

        /**
         * A variable that stores the currently used velocity file of a series.
         * It gets updated if necessary by set_current_time.
         */
        unsigned int  current_time_step;

        /**
         * Directory in which the gplates velocity are present.
         */
        std::string data_directory;

        /**
         * First part of filename of velocity files. The files have to have the pattern
         * velocity_file_name.n.gpml where n is the number of the current timestep (starts
         * from 0).
         */
        std::string velocity_file_name;

        /**
         * Time in model units (depends on other model inputs) between two velocity files.
         */
        double time_step;

        /**
         * Weight between velocity file n and n+1 while the current time is between the
         * two values t(n) and t(n+1).
         */
        double time_weight;

        /**
         * State whether we have time_dependent boundary conditions. Switched off after
         * finding no more velocity files to suppress attempts to read in new files.
         */
        bool time_dependent;

        /**
         * Pointer to an object that reads and processes data we get from gplates files.
         */
        std_cxx1x::shared_ptr<internal::GPlatesLookup> lookup;

        /**
         * Create a filename out of the name template.
         */
        std::string
        create_filename (const int timestep) const;
    };
  }
}


#endif
