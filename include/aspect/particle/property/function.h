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

#ifndef _aspect_particle_property_function_h
#define _aspect_particle_property_function_h

#include <aspect/particle/property/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * A class that initializes particle properties based on a
       * functional description provided in the input file.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class Function : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          Function();

          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;


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
           * A function object representing the particle property.
           */
          std::unique_ptr<Functions::ParsedFunction<dim>> function;

          /**
           * The coordinate representation to evaluate the function. Possible
           * choices are depth, cartesian and spherical.
           */
          Utilities::Coordinates::CoordinateSystem coordinate_system;

          /**
           * A private variable that stores the number of particle property
           * function components.
           *
           * This variable is read from the parameter file through a parameter
           * called 'Number of components'.
           */
          unsigned int n_components;
      };
    }
  }
}

#endif
