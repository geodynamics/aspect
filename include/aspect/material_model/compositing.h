/*
  Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_averaging_h
#define _aspect_material_model_averaging_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Property
    {
      /**
       * An enum to define what material properties are available.
       */
      enum  MaterialProperty
      {
        viscosity,
        density,
        thermal_expansion_coefficient,
        specific_heat,
        thermal_conductivity,
        compressibility,
        entropy_derivative_pressure,
        entropy_derivative_temperature,
        reaction_terms
      };

      static const std::map<std::string, MaterialProperty>
          property_map = {
              {"Viscosity", viscosity},
              {"Density", density},
              {"Thermal expansion coefficient", thermal_expansion_coefficient},
              {"Specific heat", specific_heat},
              {"Thermal conductivity", thermal_conductivity},
              {"Compressibility", compressibility},
              {"Entropy derivative pressure", entropy_derivative_pressure},
              {"Entropy derivative temperature", entropy_derivative_temperature},
              {"Reaction terms", reaction_terms}
          };
    }

    /**
     * A material model that selects properties from other (non-averaging) material models.
     * @ingroup MaterialModels
     */

    template <int dim>
    class Compositing : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * @copydoc MaterialModel::Interface::evaluate()
         */
        virtual
        void
        evaluate (const typename Interface<dim>::MaterialModelInputs &in,
                  typename Interface<dim>::MaterialModelOutputs &out) const;

        /**
         * @copydoc MaterialModel::Interface::declare_parameters()
         */
        static void
        declare_parameters (ParameterHandler &prm);

        /**
         * @copydoc MaterialModel::Interface::parse_parameters()
         */
        virtual void
        parse_parameters (ParameterHandler &prm);

        /**
         * @copydoc MaterialModel::Interface::is_compressible()
         *
         * Returns value from material model providing compressibility.
         */
        virtual bool is_compressible () const;

        /**
         * @copydoc MaterialModel::Interface::reference_viscosity()
         *
         * Taken from the material model providing viscosities
         */
        virtual double reference_viscosity () const;



      private:
        /**
         * Parse a string representing one of the components provided by the material models,
         * and return the resulting MaterialProperty value.
         */
        Property::MaterialProperty
        parse_property_name(const std::string &s);

        /**
         * Obtain the output properties from the component material models
         */
        void
        composite(const unsigned int model_index,
                  const typename Interface<dim>::MaterialModelOutputs &evaluated,
                  typename Interface<dim>::MaterialModelOutputs &out) const;

        /**
         * Map specifying which components a material model is responsible for
         */
        std::map<Property::MaterialProperty, unsigned int> model_property_map;

        /**
         * Pointers to the material models used for compositing
         */
        std::vector<std::string> model_names;
        std::vector<std_cxx11::shared_ptr<Interface<dim> > > models;
    };
  }
}

#endif
