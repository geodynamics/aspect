/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_ascii_data_slice_h
#define _aspect_initial_composition_ascii_data_slice_h

#include <aspect/initial_composition/interface.h>
#include <aspect/boundary_velocity/gplates.h>

#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements the prescribed compositional fields determined
     * from a AsciiData input file. Uses only a 2d slice from a 3d ascii data
     * input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class AsciiDataSlice : public Interface<dim>, public Utilities::AsciiDataInitial<dim, 3>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiDataSlice ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim, 3>::initialize;

        /**
         * Return the initial composition as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const override;

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
        // The point that defines the location of the slice.
        Tensor<1,3> slice_normal_vector;

        /**
          * The matrix, which describes the rotation by which a 2D model
          * needs to be transformed to a plane that contains the origin and
          * the prescribed point given in the input.
          */
        Tensor<2,3> rotation_matrix;

        /**
         * A function that returns the corresponding euler angles for a
         * rotation described by rotation axis and angle.
         */
        Tensor<2,3>
        rotation_matrix_from_axis (const Tensor<1,3> &rotation_axis,
                                   const double rotation_angle) const;
    };
  }
}


#endif
