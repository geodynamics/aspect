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


#ifndef _aspect_mesh_refinement_artificial_viscosity_h
#define _aspect_mesh_refinement_artificial_viscosity_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion based on the
     * artificial viscosity of the temperature or compositional fields.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class ArtificialViscosity : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * @copydoc Interface<dim>::execute()
         */
        void
        execute (Vector<float> &error_indicators) const override;

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
         * Scaling factor for the temperature indicator.
         *
         * This variable is read from the parameter file through a parameter called 'Temperature scaling factor'.
         */
        double temperature_scaling_factor;

        /**
         * The scaling factors for each compositional field.
         *
         * This variable is read from the parameter file through a parameter called 'Compositional field scaling factors'.
         */
        std::vector<double> composition_scaling_factors;
    };
  }
}

#endif
