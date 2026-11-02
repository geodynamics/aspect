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
