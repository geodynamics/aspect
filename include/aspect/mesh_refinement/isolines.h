/*
  Copyright (C) 2014 - 2019 by the authors of the ASPECT code.

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



#ifndef _aspect_isolines_h
#define _aspect_isolines_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MeshRefinement
  {
    namespace Internal
    {
      enum class PropertyName
      {
        Temperature,
        Composition
      };

      class Property
      {
        public:
          PropertyName name;
          unsigned int index;
      };

      class Isoline
      {
        public:
          /**
           * Todo
           */
          bool values_are_in_range(const std::vector<double> values) const;

          std::vector<double> min_values;
          std::vector<double> max_values;
          double min_refinement;
          double max_refinement;
          std::vector<Property> properties;
      };


    }

    /**
     * A class that implements a Isolines TODO
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Isolines : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * At the beginning of each time step, update the time for the
         * ParsedFunction.
         */
        void
        update () override;

        /**
         * After cells have been marked for coarsening/refinement, apply
         * additional criteria independent of the error estimate.
         *
         */
        void
        tag_additional_cells () const override;

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
         * Todo
         */
        std::vector<Internal::Isoline> isolines;

    };
  }
}

#endif
