/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
*/


#ifndef __aspect__initial_conditions_polygon_depth_h
#define __aspect__initial_conditions_polygon_depth_h

#include <aspect/initial_conditions/interface.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialConditions
  {
    using namespace dealii;

    /**
     * A class that describes the world initial temperature field for
     * a ellipsoidal chunk geometry.
     *
     * @ingroup InitialConditionsModels
     */
    template <int dim>
    class PolygonDepth : public InitialConditions::Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        virtual
        double initial_temperature(const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file
        */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double rad_to_degree = 180/numbers::PI;
        double degree_to_rad = numbers::PI/180;

        double reference_temperature, bottom_depth;
        std::vector<Point<2> > coordinate_list;


        /*WorldGeneratorParameters<dim> world_generator_parameters;

        double specific_heat_Cp;
        double thermal_conductivity_k;
        double thermal_expansion_coefficient_alfa;
        double thermal_diffusivity_kappa;
        double density;
        double surface_temperature;
        double potential_mantle_temperature;
        bool has_sticky_air;
        unsigned int number_of_objects;
        std::vector<int> number_of_coordinates;
        //std::vector<int> number_of_radii;
        std::vector<std::string> module_name;
        std::vector<std::string> object_name;
        std::vector<std::string> temperature_submodule;
        std::vector<std::string> composition_submodule;
        std::vector<std::string> temperature_submodule_parameters;
        std::vector<std::string> temperature_string_cache; // may give problems because it may be possible that it is shared over multiple processors causing the cache to become corrupted...
        std::vector<std::vector<std::string> > temperature_vector_cache;
        std::vector<std::vector<std::string> > temperature_submodule_vector_parameters;
        std::vector<std::vector<std::string> > coordinate_list_strings;
        std::vector<std::vector<Point<2> > > coordinate_list;
        std::vector<std::vector<std::vector<std::string> > > temperature_submodule_vector_vector_parameters;
        std::vector<std::vector<std::vector<std::string> > > temperature_2dvector_cache;
        //std::vector<std::vector<std::string> > radii_list;
        std::vector<Utilities::tk::spline> x_spline;
        std::vector<Utilities::tk::spline> y_spline;*/
    };
  }
}

#endif

