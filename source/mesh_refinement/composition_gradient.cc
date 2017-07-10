/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/composition_gradient.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    CompositionGradient<dim>::execute(Vector<float> &indicators) const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("This refinement criterion can not be used when no "
                               "compositional fields are active!"));

      indicators = 0;
      Vector<float> this_indicator (indicators.size());
      const double power = 1.0 + dim/2.0;

      const Quadrature<dim> quadrature(this->get_fe().base_element(this->introspection().base_elements.compositional_fields).get_unit_support_points());
      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_gradients);

      // the values of the compositional fields are stored as block vectors for each field
      // we have to extract them in this structure
      std::vector<Tensor<1,dim> > composition_gradients (quadrature.size());

      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();
          unsigned int i=0;
          for (; cell!=endc; ++cell, ++i)
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);
                fe_values[this->introspection().extractors.compositional_fields[c]].get_function_gradients (this->get_solution(),
                    composition_gradients);

                // For each composition dof, write into the output vector the
                // composition gradient. Note that quadrature points and dofs
                // are enumerated in the same order.
                for (unsigned int j=0; j<this->get_fe().base_element(this->introspection().base_elements.compositional_fields).dofs_per_cell; ++j)
                  this_indicator[i] += composition_gradients[j].norm();

                // Scale gradient in each cell with the correct power of h. Otherwise,
                // error indicators do not reduce when refined if there is a density
                // jump. We need at least order 1 for the error not to grow when
                // refining, so anything >1 should work. (note that the gradient
                // itself scales like 1/h, so multiplying it with any factor h^s, s>1
                // will yield convergence of the error indicators to zero as h->0)
                this_indicator[i] *= std::pow(cell->diameter(), power);
              }
          indicators.add(composition_scaling_factors[c], this_indicator);
        }
    }

    template <int dim>
    void
    CompositionGradient<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition gradient");
        {
          prm.declare_entry("Compositional field scaling factors",
                            "",
                            Patterns::List (Patterns::Double(0)),
                            "A list of scaling factors by which every individual compositional "
                            "field gradient will be multiplied. If only a single compositional "
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
    CompositionGradient<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition gradient");
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(CompositionGradient,
                                              "composition gradient",
                                              "A mesh refinement criterion that computes refinement "
                                              "indicators from the gradients of compositional fields. "
                                              "If there is more than one compositional field, then "
                                              "it simply takes the sum of the indicators times "
                                              "a user-specified weight for each field."
                                              "\n\n"
                                              "This refinement criterion computes the gradient of "
                                              "the compositional field at quadrature points on each "
                                              "cell, and then averages them in some way to obtain a "
                                              "refinement indicator for each cell. This will give a "
                                              "reasonable approximation of the true gradient of the "
                                              "compositional field if you are using a continuous "
                                              "finite element."
                                              "\n\n"
                                              "On the other hand, for discontinuous "
                                              "finite elements (see the `Use discontinuous composition "
                                              "discretization' parameter in the `Discretization' "
                                              "section), the gradient at quadrature points does not "
                                              "include the contribution of jumps in the compositional "
                                              "field between cells, and consequently will not be an "
                                              "accurate approximation of the true gradient. As an "
                                              "extreme example, consider the case of using piecewise "
                                              "constant finite elements for compositional fields; in "
                                              "that case, the gradient of the solution at quadrature "
                                              "points inside each cell will always be exactly zero, "
                                              "even if the finite element solution is different "
                                              "from each cell to the next. Consequently, the "
                                              "current refinement criterion will likely not be "
                                              "useful in this situation. That said, "
                                              "the `composition approximate "
                                              "gradient' refinement criterion exists for exactly "
                                              "this purpose.")
  }
}
