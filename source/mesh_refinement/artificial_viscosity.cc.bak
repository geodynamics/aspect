/*
  Copyright (C) 2015 - 2020 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/artificial_viscosity.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    ArtificialViscosity<dim>::execute(Vector<float> &indicators) const
    {
      indicators = 0;
      Vector<float> this_indicator(indicators.size());
      if (temperature_scaling_factor > 0.0)
        {
          this->get_artificial_viscosity(indicators);
          indicators *= temperature_scaling_factor;
        }

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          this_indicator = 0;
          this->get_artificial_viscosity_composition(this_indicator, c);

          // compute indicators += c*this_indicator:
          indicators.add(composition_scaling_factors[c], this_indicator);
        }
    }

    template <int dim>
    void
    ArtificialViscosity<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Artificial viscosity");
        {
          prm.declare_entry("Temperature scaling factor",
                            "0.0",
                            Patterns::Double (0.),
                            "A scaling factor for the artificial viscosity "
                            " of the temperature equation. Use 0.0 to disable.");
          prm.declare_entry("Compositional field scaling factors",
                            "",
                            Patterns::List (Patterns::Double (0.)),
                            "A list of scaling factors by which every individual compositional "
                            "field will be multiplied. These "
                            "factors are used to weigh the various indicators relative to "
                            "each other and to the temperature. "
                            "\n\n"
                            "If the list of scaling factors given in this parameter is empty, then this "
                            "indicates that they should all be chosen equal to 0. If the list "
                            "is not empty then it needs to have as many entries as there are "
                            "compositional fields.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    ArtificialViscosity<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Artificial viscosity");
        {
          temperature_scaling_factor = prm.get_double("Temperature scaling factor");

          composition_scaling_factors
            = Utilities::string_to_double(
                Utilities::split_string_list(prm.get("Compositional field scaling factors")));

          AssertThrow (composition_scaling_factors.size() == this->n_compositional_fields()
                       ||
                       composition_scaling_factors.size() == 0,
                       ExcMessage ("The number of scaling factors given here must either be "
                                   "zero or equal to the number of chosen refinement criteria."));

          if (composition_scaling_factors.size() == 0)
            composition_scaling_factors.resize (this->n_compositional_fields(), 0.0);

          const double sum_composition_factors =
            std::accumulate (composition_scaling_factors.begin(),
                             composition_scaling_factors.end(), 0.0);

          AssertThrow(sum_composition_factors + temperature_scaling_factor > 0.0,
                      ExcMessage("You need to have a positive scaling factor "
                                 "for at least one variable."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(ArtificialViscosity,
                                              "artificial viscosity",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators from the artificial viscosity "
                                              "of the temperature or compositional fields "
                                              "based on user specified weights.")
  }
}
