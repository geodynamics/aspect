/*
  Copyright (C) 2023 - 2024 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_slab_model_h
#define _aspect_initial_composition_slab_model_h

#include <aspect/initial_composition/interface.h>

#include <aspect/simulator_access.h>
#include <aspect/structured_data.h>


namespace aspect
{
  namespace InitialComposition
  {
    /**
     * A class that implements subducted slab geometries as a compositional
     * field determined from an input file. The file defines the depth to
     * the top of the slab and the slab thickness.
     *
     * An example model that is included is Slab2 described in
     * Hayes, G. P., Moore, G. L., Portner, D. E., Hearne, M., Flamme, H., Furtney, M.,
     * & Smoczyk, G. M. (2018). Slab2, a comprehensive subduction zone geometry model.
     * Science, 362(6410), 58-61.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class SlabModel : public Interface<dim>, public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * Return the initial composition as a function of position. For the
         * current class, this function returns 1.0 inside subducted slabs and 0.0 outside.
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

        /**
         * A function that returns the slab_boundary object to allow access
         * in other places.
         */
        const Utilities::AsciiDataBoundary<dim> &get_slab_boundary() const;

      private:
        /**
         * An object of ascii data boundary to input subducted slab depths
         * and thicknesses.
         */
        Utilities::AsciiDataBoundary<dim> slab_boundary;

        /**
         * Cache the surface boundary id to avoid unnecessary lookups.
         */
        types::boundary_id surface_boundary_id;

        /**
         * Cache the compositional field index that corresponds to
         * the slab composition to avoid unnecessary lookups.
         */
        unsigned int slab_index;
    };
  }
}


#endif
