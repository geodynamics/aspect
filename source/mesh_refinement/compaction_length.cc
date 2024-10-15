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


#include <aspect/mesh_refinement/compaction_length.h>
#include <aspect/melt.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    CompactionLength<dim>::tag_additional_cells () const
    {
      // tag_additional_cells is executed before the equations are solved
      // for the very first time. If we do not have the finite element, we
      // do not have any material properties and just do nothing in this plugin.
      if (this->get_dof_handler().n_locally_owned_dofs() == 0)
        return;

      // Use a quadrature in the support points of the porosity to compute the
      // compaction length at:
      const typename Simulator<dim>::AdvectionField porosity = Simulator<dim>::AdvectionField::composition(
                                                                 this->introspection().compositional_index_for_name("porosity"));

      const unsigned int base_element_index = porosity.base_element(this->introspection());
      const Quadrature<dim> quadrature(this->get_fe().base_element(base_element_index).get_unit_support_points());

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature,
                               update_quadrature_points | update_values | update_gradients);

      MaterialModel::MaterialModelInputs<dim> in(quadrature.size(), this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(quadrature.size(), this->n_compositional_fields());
      MeltHandler<dim>::create_material_model_outputs(out);

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            bool refine = false;
            bool clear_coarsen = false;

            fe_values.reinit(cell);
            in.reinit(fe_values, cell, this->introspection(), this->get_solution());
            this->get_material_model().evaluate(in, out);

            MaterialModel::MeltOutputs<dim> *melt_out = out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
            AssertThrow(melt_out != nullptr,
                        ExcMessage("Need MeltOutputs from the material model for computing the melt properties."));

            // for each composition dof, check if the compaction length exceeds the cell size
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                const double compaction_length = std::sqrt((out.viscosities[i] + 4./3. * melt_out->compaction_viscosities[i])
                                                           * melt_out->permeabilities[i] / melt_out->fluid_viscosities[i]);

                // If the compaction length exceeds the cell diameter anywhere in the cell, cell is marked for refinement.
                // Do not apply any refinement if the porosity is so small that melt can not migrate.
                if (compaction_length < 2.0 * cells_per_compaction_length * cell->minimum_vertex_distance()
                    && this->get_melt_handler().is_melt_cell(cell))
                  clear_coarsen = true;

                if (compaction_length < cells_per_compaction_length * cell->minimum_vertex_distance()
                    && this->get_melt_handler().is_melt_cell(cell))
                  {
                    refine = true;
                    break;
                  }
              }

            if (clear_coarsen)
              cell->clear_coarsen_flag ();
            if (refine)
              cell->set_refine_flag ();
          }
    }

    template <int dim>
    void
    CompactionLength<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Compaction length");
        {
          prm.declare_entry("Mesh cells per compaction length", "1.0",
                            Patterns::Double (0.),
                            "The desired ratio between compaction length and size of the "
                            "mesh cells, or, in other words, how many cells the mesh should "
                            "(at least) have per compaction length. Every cell where this "
                            "ratio is smaller than the value specified by this parameter "
                            "(in places with fewer mesh cells per compaction length) is "
                            "marked for refinement.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    CompactionLength<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("Compaction length");
        {
          cells_per_compaction_length = prm.get_double ("Mesh cells per compaction length");
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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(CompactionLength,
                                              "compaction length",
                                              "A mesh refinement criterion for models with melt transport "
                                              "that computes refinement indicators based on the compaction "
                                              "length, defined as "
                                              "$\\delta = \\sqrt{\\frac{(\\xi + 4 \\eta/3) k}{\\eta_f}}$. "
                                              "$\\xi$ is the bulk viscosity, $\\eta$ is the shear viscosity, "
                                              "$k$ is the permeability and $\\eta_f$ is the melt viscosity. "
                                              "If the cell width or height exceeds a multiple (which is "
                                              "specified as an input parameter) of this compaction length, "
                                              "the cell is marked for refinement.")
  }
}
