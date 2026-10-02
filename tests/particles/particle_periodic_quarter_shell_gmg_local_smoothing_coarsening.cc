/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT and is distributed under the GNU General
  Public License, version 2 or later. See the file LICENSE.
*/

#include <aspect/mesh_refinement/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/particle/manager.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  template <int dim>
  class PeriodicRefinementCycle : public MeshRefinement::Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:
      void tag_additional_cells() const override
      {
        if (this->get_dof_handler().n_dofs() == 0)
          {
            for (const auto &cell : this->get_triangulation().active_cell_iterators())
              if (cell->is_locally_owned() && cell->center()[0] > cell->center()[1])
                cell->clear_refine_flag();
            return;
          }

        const unsigned int phase = this->get_timestep_number() % 4;
        for (const auto &cell : this->get_triangulation().active_cell_iterators())
          if (cell->is_locally_owned())
            {
              cell->clear_refine_flag();
              cell->clear_coarsen_flag();
              const auto point = cell->center();
              if (phase == 0 && point[1] > std::abs(point[0]) && point.norm() > 0.75)
                cell->set_refine_flag();
              else if (phase == 2 && point[0] > std::abs(point[1]) && point.norm() < 0.75)
                cell->set_refine_flag();
              else if (phase == 1 || phase == 3)
                cell->set_coarsen_flag();
            }
      }
  };

  template <int dim>
  class PeriodicRefinementCheck : public Postprocess::Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:
      std::pair<std::string, std::string> execute(TableHandler &) override
      {
        for (const auto &cell : this->get_triangulation().active_cell_iterators())
          if (cell->is_locally_owned())
            for (const auto face : cell->face_indices())
              if (cell->has_periodic_neighbor(face))
                AssertThrow(cell->periodic_neighbor(face)->is_active()
                            && cell->periodic_neighbor(face)->level() == cell->level(),
                            ExcMessage("Periodic refinement levels differ."));

        AssertThrow(this->get_particle_manager(0).n_global_particles() == 10,
                    ExcMessage("Particles were lost during mesh adaptation."));
        const auto cells = this->get_triangulation().n_global_active_cells();
        coarsened |= cells < previous_cells;
        previous_cells = cells;
        if (this->get_timestep_number() == 4)
          AssertThrow(coarsened, ExcMessage("The test did not coarsen the mesh."));
        return {"Periodic refinement and particle checks:", coarsened ? "passed; coarsened" : "passed"};
      }

    private:
      types::global_cell_index previous_cells = 0;
      bool coarsened = false;
  };

  ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(PeriodicRefinementCycle,
                                            "periodic refinement cycle",
                                            "Alternate asymmetric refinement and coarsening.")
  ASPECT_REGISTER_POSTPROCESSOR(PeriodicRefinementCheck,
                                "periodic refinement check",
                                "Check periodic partners and particle conservation after adaptation.")
}
