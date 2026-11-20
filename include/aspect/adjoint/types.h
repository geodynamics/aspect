/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#ifndef _aspect_adjoint_types_h
#define _aspect_adjoint_types_h

#include <aspect/global.h>

#include <string>

namespace aspect
{
  namespace Adjoint
  {
    /**
     * An enum whose members identify physical material properties for which
     * adjoint kernels can currently be accumulated.
     *
     * @ingroup Adjoint
     */
    enum class PhysicalProperty
    {
      /**
       * Density.
       */
      density,

      /**
       * Viscosity.
       */
      viscosity
    };

    /**
     * Return the parameter-file name of a physical property.
     *
     * @param property The physical property.
     *
     * @return The name of @p property.
     *
     * @ingroup Adjoint
     */
    inline std::string
    property_name(const PhysicalProperty property)
    {
      switch (property)
        {
          case PhysicalProperty::density:
            return "density";
          case PhysicalProperty::viscosity:
            return "viscosity";
          default:
            AssertThrow(false, ExcMessage("Unknown adjoint physical property."));
            return "unknown";
        }
    }

    /**
     * A key that identifies one adjoint kernel contribution by objective,
     * physics term, and physical property.
     *
     * @ingroup Adjoint
     */
    struct KernelContributionKey
    {
      /**
       * Objective functional name associated with the contribution.
       */
      std::string objective_name;

      /**
       * Name of the physics term associated with the contribution.
       */
      std::string physics_term_name;

      /**
       * Physical property associated with the contribution.
       */
      PhysicalProperty property;

      /**
       * Compare contribution keys for use in associative containers.
       *
       * @param other The key to compare against.
       *
       * @return true if this key should be ordered before @p other.
       */
      bool
      operator<(const KernelContributionKey &other) const
      {
        if (objective_name != other.objective_name)
          return objective_name < other.objective_name;

        if (physics_term_name != other.physics_term_name)
          return physics_term_name < other.physics_term_name;

        return static_cast<unsigned int>(property) < static_cast<unsigned int>(other.property);
      }
    };
  }
}

#endif
