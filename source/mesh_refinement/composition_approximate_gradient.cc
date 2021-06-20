/*
  Copyright (C) 2015 - 2021 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/composition_approximate_gradient.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/derivative_approximation.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    CompositionApproximateGradient<dim>::execute(Vector<float> &indicators) const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("This refinement criterion can not be used when no "
                               "compositional fields are active!"));
      indicators = 0;
      // create a vector with the requisite ghost elements
      // and use it for estimating the gradients of compositional field
      LinearAlgebra::BlockVector vec(this->introspection().index_sets.system_partitioning,
                                     this->introspection().index_sets.system_relevant_partitioning,
                                     this->get_mpi_communicator());
      Vector<float> indicators_tmp(this->get_triangulation().n_active_cells());
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        {
          const unsigned int block_idx = this->introspection().block_indices.compositional_fields[c];
          vec.block(block_idx) = this->get_solution().block(block_idx);
          vec.compress(VectorOperation::insert);
          indicators_tmp = 0;
          DerivativeApproximation::approximate_gradient(this->get_mapping(),
                                                        this->get_dof_handler(),
                                                        vec,
                                                        indicators_tmp,
                                                        this->introspection().component_indices.compositional_fields[c]);

          indicators_tmp *= composition_scaling_factors[c];

          // Scale approximated gradient in each cell with the correct power of h. Otherwise,
          // error indicators do not reduce when refined if there is a density
          // jump. We need at least order 1 for the error not to grow when
          // refining, so anything >1 should work. (note that the gradient
          // itself scales like 1/h, so multiplying it with any factor h^s, s>1
          // will yield convergence of the error indicators to zero as h->0)
          const double power = 1.0 + dim / 2.0;
          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              indicators_tmp(cell->active_cell_index()) *= std::pow(cell->diameter(), power);

          indicators += indicators_tmp;
        }
    }


    template <int dim>
    void
    CompositionApproximateGradient<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition approximate gradient");
        {
          prm.declare_entry("Compositional field scaling factors",
                            "",
                            Patterns::List (Patterns::Double (0.)),
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
    CompositionApproximateGradient<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition approximate gradient");
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(CompositionApproximateGradient,
                                              "composition approximate gradient",
                                              "A mesh refinement criterion that computes refinement "
                                              "indicators from the gradients of compositional fields. "
                                              "If there is more than one compositional field, then "
                                              "it simply takes the sum of the indicators times "
                                              "a user-specified weight for each field."
                                              "\n\n"
                                              "In contrast to the `composition gradient' refinement "
                                              "criterion, the current criterion does not "
                                              "compute the gradient at quadrature points on each cell, "
                                              "but by a finite difference approximation between the "
                                              "centers of cells. Consequently, it also works if "
                                              "the compositional fields are computed using discontinuous "
                                              "finite elements.")
  }
}
