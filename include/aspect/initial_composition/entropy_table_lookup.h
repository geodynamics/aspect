// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2017-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


#ifndef _aspect_initial_composition_entropy_table_lookup_h
#define _aspect_initial_composition_entropy_table_lookup_h

#include <aspect/initial_composition/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements initial conditions for the entropy field
     * Note that this plugin only
     * works if there is a compositional field called 'entropy'.
     * All compositional fields except entropy are not changed by this plugin.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class EntropyTableLookUp : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the plugin.
         */
        void initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const override;

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
         * Location of data files. This variable is read from the parameter
         * file through the 'Data directory' parameter.
         */
        std::string data_directory;
        /**
         * This variable is read from the parameter file through a parameter called 'Material file name'.
         */
        std::string material_file_name;

        /**
         * Index of the entropy in the compositional fields
         */
        unsigned entropy_index;

        /**
         * A shared pointer to the initial temperature object
         * that ensures that the current object can continue
         * to access the initial temperature object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialTemperature::Manager<dim>> initial_temperature_manager;

        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition_manager;

        /**
         * Pointer to the StructuredDataLookup object that holds the material data.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<2>> material_lookup;
    };
  }
}


#endif
