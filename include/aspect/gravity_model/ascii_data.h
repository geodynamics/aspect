/*
  Copyright (C) 2016 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_gravity_model_ascii_data_h
#define _aspect_gravity_model_ascii_data_h

#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace GravityModel
  {
    /**
     * A class that implements a gravity description based on
     * an AsciiData input file.
     *
     * @ingroup GravityModels
     */
    template <int dim>
    class AsciiData : public Utilities::AsciiDataProfile<dim>, public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        AsciiData ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        // avoid -Woverloaded-virtual:
        using Utilities::AsciiDataProfile<dim>::initialize;

        /**
         * Return the gravity as a function of position.
         */
        Tensor<1,dim> gravity_vector (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Object containing the data profile.
         */
        std::unique_ptr<aspect::Utilities::StructuredDataLookup<1>> profile;

        /**
         * The column index of the gravity column in the data file.
         */
        unsigned int gravity_index;
    };
  }
}


#endif
