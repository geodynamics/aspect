/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/composition_threshold.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    CompositionThreshold<dim>::tag_additional_cells () const
    {
      AssertThrow (this->n_compositional_fields() >= 1,
                   ExcMessage ("This refinement criterion can not be used when no "
                               "compositional fields are active!"));

      // tag_additional_cells is executed before the equations are solved
      // for the very first time. If we do not have the finite element, we
      // do not have the compositional fields and just do nothing in this plugin.
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;

      const unsigned int dofs_per_cell = this->get_fe().dofs_per_cell;
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            cell->get_dof_indices (local_dof_indices);
            bool refine = false;

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              {
                const unsigned int component_idx = this->introspection().component_indices.compositional_fields[c];
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  {
                    if (this->get_fe().system_to_component_index(i).first == component_idx)
                      {
                        const double composition_value = this->get_solution()[local_dof_indices[i]];
                        // if the composition exceeds the threshold, cell is marked for refinement
                        if (composition_value > composition_thresholds[c])
                          {
                            refine = true;
                            break;
                          }
                      }
                  }
                if (refine)
                  break;
              }

            if (refine)
              {
                cell->clear_coarsen_flag ();
                cell->set_refine_flag ();
              }
          }
    }

    template <int dim>
    void
    CompositionThreshold<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition threshold");
        {
          prm.declare_entry("Compositional field thresholds",
                            "",
                            Patterns::List (Patterns::Double()),
                            "A list of thresholds, one for each compositional field "
                            "to be evaluated against.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    CompositionThreshold<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Composition threshold");
        {
          composition_thresholds
            = Utilities::string_to_double(
                Utilities::split_string_list(prm.get("Compositional field thresholds")));

          AssertThrow (composition_thresholds.size() == this->n_compositional_fields(),
                       ExcMessage ("The number of thresholds given here must be "
                                   "equal to the number of compositional fields."));
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(CompositionThreshold,
                                              "composition threshold",
                                              "A mesh refinement criterion that computes refinement "
                                              "indicators from the compositional fields. One threshold "
                                              "per compositional is given in the input file, and if any "
                                              "field exceeds its threshold, the cell is marked for "
                                              "refinement.")
  }
}
