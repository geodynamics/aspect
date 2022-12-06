/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_ascii_data_h
#define _aspect_initial_composition_ascii_data_h

#include <aspect/initial_composition/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements the prescribed compositional fields determined
     * from a AsciiData input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class AsciiData : public Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Empty Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

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
        /**
         * The dataset that contains the initial composition.
         * This member variable is used in case the dataset has the same
         * spatial dimensions as the model.
         */
        std::unique_ptr<Utilities::AsciiDataInitial<dim>> ascii_data_initial;

        /**
         * The dataset that contains the initial composition.
         * This member variable is used in case the dataset has a higher
         * spatial dimension than the model, and only a slice
         * of the dataset is used.
         */
        std::unique_ptr<Utilities::AsciiDataInitial<dim,3>> ascii_data_slice;

        /**
         * Whether to use a dataset that has the same spatial dimensions as
         * the model or not. If true only a slice of the dataset is used.
         */
        bool slice_data;

        /**
         * The two points that define the location of the slice.
         */
        Tensor<1,3> first_point_on_slice;
        Tensor<1,3> second_point_on_slice;

        /**
         * The matrix that describes the rotation by which a 2D model
         * needs to be transformed to a plane that contains the origin and
         * the two prescribed points given in the input.
         */
        Tensor<2,3> rotation_matrix;
    };
  }
}


#endif
