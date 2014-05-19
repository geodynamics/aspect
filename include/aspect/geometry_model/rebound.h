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


#ifndef __aspect__geometry_model_rebound_h
#define __aspect__geometry_model_rebound_h

#include <aspect/geometry_model/interface.h>
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace GeometryModel
  {
    using namespace dealii;

    /**
     * A class deriving from Box<dim>, which changes the upper boundary 
       with a sinusoidal perturbation of a given order and amplitude
     */
    template <int dim>
    class ReboundBox : public Box<dim>
    {
      public:
        /**
         * Generate a coarse mesh for the geometry described by this class.
         * Makes perturbs the top boundary of the box with a function
         * of the form z' = amplitude * cos(order * x )
         */
        virtual
        void create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const;

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
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:

        unsigned int order;  //Order of the perturbation
        double amplitude;  //amplitude of the perturbation

    };
  }
}

#endif
