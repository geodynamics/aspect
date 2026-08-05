/*
  Copyright (C) 2026 by the authors of the ASPECT code.

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


#include <aspect/postprocess/finite_element_information.h>

#include <limits>
#include <sstream>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    FiniteElementInformation<dim>::execute (TableHandler &)
    {
      // Only report once:
      if (this->get_timestep_number() > 0 || this->get_pre_refinement_step() != std::numeric_limits<unsigned int>::max())
        return {"", ""};

      std::ostringstream output;

      unsigned int compositional_field_index = 0;
      const auto &compositional_field_names = this->introspection().get_composition_names();
      const auto &compositional_field_descriptions = this->introspection().get_composition_descriptions();

      for (const auto &variable : this->introspection().get_variables())
        if (variable.name == "compositions")
          for (unsigned int i = 0; i < variable.multiplicity; ++i)
            {
              output << "composition " << compositional_field_index << ": "
                     << compositional_field_names[compositional_field_index]
                     << " (" << CompositionalFieldDescription::type_to_string(compositional_field_descriptions[compositional_field_index].type)
                     << "): " << variable.fe->get_name()
                     << std::endl;
              ++compositional_field_index;
            }
        else
          output << variable.name << ": "
                 << variable.fe->get_name() << std::endl;

      return {"Finite element spaces:", output.str()};
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(FiniteElementInformation,
                                  "finite element information",
                                  "A postprocessor that prints the names and finite element "
                                  "spaces of all solution variables at time zero. "
                                  "For compositional fields, it also prints the compositional "
                                  "field type.")
  }
}
