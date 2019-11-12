/*
  Copyright (C) 2016 - 2019-2018 by the authors of the ASPECT code.

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


#ifndef _aspect_mesh_refinement_volume_of_fluid_interface_h
#define _aspect_mesh_refinement_volume_of_fluid_interface_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion that refines the
     * mesh near boundaries for the volume of fluid interface tracking algorithm.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class VolumeOfFluidInterface: public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Mark large unrefined neighboring cells for refinement and prevent
         * coarsening
         */
        void
        tag_additional_cells() const override;

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

      private:
        /**
         * If true, flag cells for coarsening if they are not flagged for refinement
         */
        bool strict_coarsening;
    };
  }
}

#endif
