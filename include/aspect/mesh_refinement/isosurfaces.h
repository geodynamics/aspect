/*
  Copyright (C) 2014 - 2021 by the authors of the ASPECT code.

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



#ifndef _aspect_mesh_refinement_isosurfaces_h
#define _aspect_mesh_refinement_isosurfaces_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace MeshRefinement
  {
    namespace internal
    {
      enum class PropertyType
      {
        Temperature,
        Composition
      };

      struct Property
      {
        public:
          /**
           * Constructor. Converts a property name into a structure containing a property type
           * and an index. If the property contains multiple items (e.g. the property compositional field has a
           * field index) the index referring to the particular item of that property is stored in the variable index.
           * @param property_name The name of a property, which can be Temperature for the temperature field or
           * the name of a compositional field listed in the parameter available_compositions.
           * @param available_compositions A list of names of the available compositional fields.
           */
          Property(const std::string &property_name,
                   const std::vector<std::string> &available_compositions);

          /**
           * The Property type of the property
           */
          PropertyType type;

          /**
           * An index, in case the property type contains multiple values. This is
           * currently only used for storing which compositional field the property
           * corresponds to.
           */
          unsigned int index;
      };

      struct Isosurface
      {
        public:
          /**
           * Checks whether the provided values are between the min and max value for the
           * set properties of the isosurface. This function assumes that the order of the
           * provided @p values matches the order in which the properties are stored.
           *
           * This function assumes that @p values and all the vectors in isosurfaces already
           * have the same length.
           */
          bool are_all_values_in_range(const std::vector<double> &values) const;

          std::vector<double> min_values;
          std::vector<double> max_values;
          int min_refinement;
          int max_refinement;
          std::vector<Property> properties;
      };

    }

    /**
     * A class that implements an Isosurfaces mesh refinement plugin. This
     * plugin allows for setting a minimum and a maximum refinement level in
     * a part of the model domain where a variable/property (e.g. Temperature)
     * is between two values (e.g. two isotherms of 274K and 1600K). This is
     * currently implemented for temperature and compositions.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class Isosurfaces : public Interface<dim>,
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
         * Declare the parameters this class takes from input files.
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
         * A vector of the isosurfaces used by this class.
         */
        std::vector<internal::Isosurface> isosurfaces;

    };
  }
}

#endif
