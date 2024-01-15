/*
  Copyright (C) 2016 by the authors of the ASPECT code.

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


#ifndef _aspect_geometry_model__initial_topography_lithosphere_rift_h
#define _aspect_geometry_model__initial_topography_lithosphere_rift_h

#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace InitialTopographyModel
  {
    using namespace dealii;

    /**
     * A class that describes an initial topography for the geometry model,
     * by defining a set of polylines on the surface from the prm file. It
     * sets the elevation in each Polyline to a constant value.
     */
    template <int dim>
    class LithosphereRift : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        /**
         * Return the value of the topography for a point.
         */
        virtual
        double value (const Point<dim-1> &p) const;

        /**
         * Return the maximum value of the elevation.
         */
        virtual
        double max_topography () const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * The standard deviation of the Gaussian amplitude of the lithospheric thicknesses
         * with distance from the rift axis.
         */
        double sigma_rift;

        /**
         * The half width of the hyperbolic tangent used to smooth the transitions
         * between reference and polygon lithospheric thicknesses.
         */
        double sigma_polygon;

        /**
         * The maximum amplitude of the Gaussian distribution of the lithospheric thicknesses
         * with distance from the rift axis. It should have values between -1 and 1, where positive
         * numbers represent a reduction in thickness and negative numbers an increase.
         */
        std::vector<double> A;

        /**
         * The maximum amplitude of the Gaussian distribution of the topography around the rift.
         */
        double topo_rift_amplitude;

        /**
         * The reference lithospheric column used in computing the topography based on isostasy
         * and the thickness of this column.
         */
        double ref_rgh;
        double compensation_depth;

        /**
         * Wheter or not to take the polygon thicknesses as dominant, or to smooth them
         * gradually into rift areas.
         */
        bool blend_rift_and_polygon;


        /**
         * The maximum amplitude of the topography of the polygon area.
         */
        double topo_polygon_amplitude;

        /**
         * Vector for field densities.
         */
        std::vector<double> densities;
        std::vector<double> temp_densities;

        /**
         * Vector for the reference field thicknesses away from the rift.
         */
        std::vector<double> thicknesses;

        /**
         * Vector for the rift  thicknesses at the center.
         */
        std::vector<double> rift_thicknesses;

        /**
         * Vector for the polygon thicknesses at the center.
         */
        std::vector<std::vector<double>> polygon_thicknesses;

        /**
         * The maximum topography in this model
         */
        double maximum_topography;

        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition_manager;

        /**
         * Vector containing the number of phases for each composition.
         */
        std::unique_ptr<std::vector<unsigned int>> n_phases_for_each_composition;

        /**
         * Vector containing the names of the compositional fields.
         */
        std::vector<std::string> list_of_composition_names;
    };
  }
}


#endif
