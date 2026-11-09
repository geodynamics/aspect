// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2016-2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_mesh_deformation_ascii_data_h
#define _aspect_mesh_deformation_ascii_data_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * A class that implements initial topography determined
     * from an AsciiData input file.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataBoundary<dim>, public Interface<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Sets the boundary id of the surface boundary.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataBoundary<dim>::initialize;

        /**
         * Return the surface topography as a function of position along the surface.
         * For the current class, this function returns a value from the text files.
         *
         * @copydoc aspect::InitialTopographyModel::Interface::value()
         */
        Tensor<1,dim>
        compute_initial_deformation_on_boundary(const types::boundary_id boundary_indicator,
                                                const Point<dim> &position) const override;

        /**
         * Returns whether or not the plugin requires surface stabilization
         */
        bool needs_surface_stabilization () const override;

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
        types::boundary_id surface_boundary_id;
    };
  }
}


#endif
