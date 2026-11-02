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


#ifndef _aspect_geometry_model_initial_topography_model_ascii_data_h
#define _aspect_geometry_model_initial_topography_model_ascii_data_h

#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
     * A class that implements topography determined
     * from an AsciiData input file.
     *
     * @ingroup InitialTopographyModel
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
        double
        value (const Point<dim-1> &surface_point) const override;

        /**
         * Return the maximum value of the elevation.
         */
        double max_topography () const override;

        /**
         * Return the gradient of the surface topography for a given position
         * along the surface.
         */
        Tensor<1,dim-1>
        vector_gradient(const Point<dim> &p) const;

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
