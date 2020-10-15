/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef _aspect_mesh_deformation_ascii_data_h
#define _aspect_mesh_deformation_ascii_data_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MeshDeformation
  {
    using namespace dealii;

    /**
     * A class that implements initial topography determined
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
        Tensor<1,dim>
        compute_initial_deformation_on_boundary(const types::boundary_id boundary_indicator,
                                                const Point<dim> &position) const override;

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
