// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------


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
