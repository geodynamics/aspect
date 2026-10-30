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
         *
         * This variable is read from the parameter file through a parameter called 'Strict coarsening'.
         */
        bool strict_coarsening;
    };
  }
}

#endif
