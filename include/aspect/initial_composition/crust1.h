/*
  Copyright (C) 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_composition_crust1_h
#define _aspect_initial_composition_crust1_h

#include <aspect/initial_composition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements an initial compositional field containing
     * densities from the CRUST 1.0 data set. Note that this plugin only
     * works if there is a compositional field called 'crust1_densities'.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Crust1 : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const;


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
         * File directory and names
         */
        std::string data_directory;
        std::string crust1_bnds_file_name;
        std::string crust1_rho_file_name;

        std::array<std::array<double, 9>, 64800> bnds;
        std::array<std::array<double, 9>, 64800> rho;

        std::array<std::array<double, 2>, 64800> crust1_lon_lat;
        
    };
  }
}


#endif
