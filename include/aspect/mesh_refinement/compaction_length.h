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


#ifndef _aspect_mesh_refinement_compaction_length_h
#define _aspect_mesh_refinement_compaction_length_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion based on the
     * compaction length, a typical length scale in models with melt
     * migration. The user can specify a desired ratio between the
     * compaction length and the size of the mesh cells, and every cell
     * where this ratio is smaller than specified is marked for refinement.
     * This means that there will always be a given (minimum) number of
     * mesh cells per compaction length.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class CompactionLength : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         */
        void
        tag_additional_cells () const override;

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
         * The desired ratio between compaction length and size of the mesh cells,
         * in other words, how many cells the mesh should have per compaction length.
         *
         * This variable is read from the parameter file through a parameter called 'Mesh cells per compaction length'.
         */
        double cells_per_compaction_length;
    };
  }
}

#endif
