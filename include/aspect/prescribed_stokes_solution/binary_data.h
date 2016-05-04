/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#ifndef __aspect__perscribed_stokes_solution_binary_data_h
#define __aspect__perscribed_stokes_solution_binary_data_h

#include <deal.II/distributed/tria.h>

#include <aspect/prescribed_stokes_solution/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace PrescribedStokesSolution
  {
    using namespace dealii;

    /**
     * A class that implements a prescribed velocity field determined from
     * a BinaryData input file.
     *
     * @ingroup PrescribedStokesSolution
     */
    template <int dim>
    class BinaryData : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        BinaryData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        virtual
        void
        update();

        virtual
        void
        stokes_solution (const Point<dim> &position, Vector<double> &value) const;

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
        parse_parameters (ParameterHandler &prm);
      //private:
        std::string dataDirectory;
        //parallel::distrtibuted::Triangulation<dim> triangulation;
        parallel::distributed::Triangulation<dim> *triangulation;
        //DoFHandler<dim> dof_handler;
        LinearAlgebra::BlockVector solution;
    };
  }
}


#endif
