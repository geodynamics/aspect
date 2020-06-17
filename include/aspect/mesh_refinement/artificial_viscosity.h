/*
  Copyright (C) 2015 - 2019 by the authors of the ASPECT code.

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
         */
        double temperature_scaling_factor;

        /**
         * The scaling factors for each compositional field.
         */
        std::vector<double> composition_scaling_factors;
    };
  }
}

#endif
