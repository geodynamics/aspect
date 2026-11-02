// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2025 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_postprocess_current_surface_h
#define _aspect_postprocess_current_surface_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/base/function_lib.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that creates a function to store
     * the current surface including mesh deformation.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class CurrentSurface : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Function that uses the stored surface topography to calculate the depth
         * including changes from mesh deformation. This will linearly interpolate
         * between the two nearest surface points to get a depth for the given
         * position.
         */
        double depth_including_mesh_deformation(const Point<dim> &position) const;

        /**
         * Evaluate the solution and compute the current surface elevations.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

      private:
        /**
         * Store the coordinates for the surface function.
         */
        std::array<std::vector<double>, dim-1> coordinates;

        /**
         * Store a data table of the elevations.
         */
        Table<dim-1, double> data_table;

        /**
         * Function to store and interpolate surface topography.
         */
        std::unique_ptr<Functions::InterpolatedTensorProductGridData<dim-1>> surface_function;
    };
  }
}


#endif
