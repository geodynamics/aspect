/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
                   ExcMessage ("This refinement criterion can not be used when no "
                               "compositional fields are active!"));
      indicators = 0;

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          Vector<float> this_indicator (indicators.size());

          QGauss<dim-1> quadrature (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);

          KellyErrorEstimator<dim>::estimate (this->get_dof_handler(),
                                              quadrature,
                                              typename FunctionMap<dim>::type(),
                                              this->get_solution(),
                                              this_indicator,
                                              this->introspection().component_masks.compositional_fields[c],
                                              0,
                                              0,
                                              this->get_triangulation().locally_owned_subdomain());
          for (unsigned int i=0; i<indicators.size(); ++i)
            this_indicator[i] *= composition_scaling_factors[c];
          indicators += this_indicator;
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
                            Patterns::List (Patterns::Double(0)),
                            "A list of scaling factors by which every individual compositional "
                            "field will be multiplied by. If only a single compositional "
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
                                              "from each of the compositional field.")
  }
}
