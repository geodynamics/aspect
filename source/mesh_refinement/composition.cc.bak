/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/composition.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    Composition<dim>::execute(Vector<float> &indicators) const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("This refinement criterion cannot be used when no "
                               "compositional fields are active!"));
      indicators = 0;
      Vector<float> this_indicator (indicators.size());

      const Quadrature<dim-1> &quadrature = this->introspection().face_quadratures.compositional_fields;

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          this_indicator = 0;
          KellyErrorEstimator<dim>::estimate (this->get_mapping(),
                                              this->get_dof_handler(),
                                              quadrature,
                                              std::map<types::boundary_id,const Function<dim>*>(),
                                              this->get_solution(),
                                              this_indicator,
                                              this->introspection().component_masks.compositional_fields[c],
                                              nullptr,
                                              0,
                                              this->get_triangulation().locally_owned_subdomain());
          // compute indicators += c*this_indicator:
          indicators.add(composition_scaling_factors[c], this_indicator);
        }
    }

    template <int dim>
    void
    Composition<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition");
        {
          prm.declare_entry("Compositional field scaling factors",
                            "",
                            Patterns::List (Patterns::Double (0.)),
                            "A list of scaling factors by which every individual compositional "
                            "field will be multiplied. If only a single compositional "
                            "field exists, then this parameter has no particular meaning. "
                            "On the other hand, if multiple criteria are chosen, then these "
                            "factors are used to weigh the various indicators relative to "
                            "each other. "
                            "\n\n"
                            "If the list of scaling factors given in this parameter is empty, then this "
                            "indicates that they should all be chosen equal to one. If the list "
                            "is not empty then it needs to have as many entries as there are "
                            "compositional fields.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    Composition<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition");
        {
          composition_scaling_factors
            = Utilities::string_to_double(
                Utilities::split_string_list(prm.get("Compositional field scaling factors")));

          AssertThrow (composition_scaling_factors.size() == this->n_compositional_fields()
                       ||
                       composition_scaling_factors.size() == 0,
                       ExcMessage ("The number of scaling factors given here must either be "
                                   "zero or equal to the number of chosen refinement criteria."));

          if (composition_scaling_factors.size() == 0)
            composition_scaling_factors.resize (this->n_compositional_fields(), 1.0);
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(Composition,
                                              "composition",
                                              "A mesh refinement criterion that computes "
                                              "refinement indicators from the compositional fields. "
                                              "If there is more than one compositional field, then "
                                              "it simply takes the sum of the indicators computed "
                                              "from each of the compositional field."
                                              "\n\n"
                                              "The way these indicators are computed is by "
                                              "evaluating the `Kelly error indicator' on each "
                                              "compositional field. This error indicator takes the "
                                              "finite element approximation of the compositional "
                                              "field and uses it to compute an approximation "
                                              "of the second derivatives of the composition for "
                                              "each cell. This approximation is then multiplied "
                                              "by an appropriate power of the cell's diameter "
                                              "to yield an indicator for how large the error "
                                              "is likely going to be on this cell. This "
                                              "construction rests on the observation that for "
                                              "many partial differential equations, the error "
                                              "on each cell is proportional to some power of "
                                              "the cell's diameter times the second derivatives "
                                              "of the solution on that cell."
                                              "\n\n"
                                              "For complex equations such as those we solve "
                                              "here, this observation may not be strictly "
                                              "true in the mathematical sense, but it often "
                                              "yields meshes that are surprisingly good.")
  }
}
