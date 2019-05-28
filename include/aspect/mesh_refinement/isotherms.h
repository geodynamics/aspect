/*
  Copyright (C) 2011-2018 by the authors of the ASPECT code.

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


#ifndef __aspect__mesh_refinement_isotherms_h
#define __aspect__mesh_refinement_isotherms_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion based on
     * isotherms.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Isotherms : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         */
        virtual
        void
        tag_additional_cells () const;

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

      private:
        /**
         * Declare a variable which indicates if and which compositional field
         * should be excluded from the process.
         */
        int exclude_composition;

        /**
         * Declare a variable which indicates how many isotherms the user has
         * declared, and a variable to store all those isotherms
         */
        std::vector<std::pair<double,double> > isotherms;

        /**
         * This variable stores the min and max refinement levels associated
         * with the isotherms variable temperatgures.
         */
        std::vector<std::pair<int,int> > isotherms_levels;

        /**
         * Declare a variable to store the global minimum refinement level and
         * a variable to store the global maximum refinement level set by the user
         */
        unsigned int minimum_refinement_level;
        unsigned int maximum_refinement_level;

    };
  }
}

#endif
