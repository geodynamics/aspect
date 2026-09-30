/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT and is distributed under the GNU General
  Public License, version 2 or later. See the file LICENSE.
*/

#include "../../benchmarks/annulus/plugin/annulus.cc"
#include <aspect/simulator/solver/stokes_matrix_free_local_smoothing.h>

#include <fstream>
#include <iomanip>

namespace aspect
{
  template <int dim>
  class PeriodicMultigridCheck : public Postprocess::Interface<dim>,
    public SimulatorAccess<dim>
  {
    public:
      std::pair<std::string, std::string> execute(TableHandler &) override
      {
        if constexpr (dim != 2)
          AssertThrow(false, ExcNotImplemented());

        const auto &solver = dynamic_cast<const StokesMatrixFreeHandlerLocalSmoothingImplementation<dim,2> &>(this->get_stokes_matrix_free());
        const auto &dofs = solver.get_velocity_multigrid_dof_handler();
        const auto &mg = solver.get_velocity_multigrid_constraints();
        const auto &transfer = solver.get_velocity_multigrid_transfer();
        const auto &fe = dofs.get_fe();
        using LevelVector = dealii::LinearAlgebra::distributed::Vector<GMGNumberType>;
        const unsigned int levels = this->get_triangulation().n_global_levels();
        std::vector<LevelVector> exact(levels);
        std::vector<AffineConstraints<double>> constraints(levels);
        double constraint_error = 0.;
        double transfer_error = 0.;

        for (unsigned int level = 0; level < levels; ++level)
          {
            const auto owned = dofs.locally_owned_mg_dofs(level);
            const auto relevant = DoFTools::extract_locally_relevant_level_dofs(dofs, level);
            exact[level].reinit(owned, relevant, this->get_mpi_communicator());
            constraints[level].reinit(owned, relevant);
            constraints[level].merge(mg.get_level_constraints(level), AffineConstraints<double>::no_conflicts_allowed, true);
            constraints[level].merge(mg.get_user_constraint_matrix(level), AffineConstraints<double>::no_conflicts_allowed, true);
            constraints[level].close();
            AssertThrow(Utilities::MPI::sum(constraints[level].n_constraints(), this->get_mpi_communicator()) > 0,
                        ExcMessage("The velocity level has no periodic constraints."));

            std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
            for (const auto &cell : dofs.mg_cell_iterators_on_level(level))
              if (!cell->is_artificial_on_level())
                {
                  cell->get_mg_dof_indices(indices);
                  for (unsigned int i = 0; i < indices.size(); ++i)
                    if (owned.is_element(indices[i]))
                      {
                        const auto point = this->get_mapping().transform_unit_to_real_cell(cell, fe.get_unit_support_points()[i]);
                        const unsigned int component = fe.system_to_component_index(i).first;
                        exact[level][indices[i]] = component == 0 ? -point[1] : point[0];
                      }
                }
            exact[level].compress(VectorOperation::insert);
            exact[level].update_ghost_values();
            const LevelVector &values = exact[level];
            for (const auto &line : constraints[level].get_lines())
              if (owned.is_element(line.index))
                {
                  double residual = values[line.index] - line.inhomogeneity;
                  for (const auto &entry : line.entries)
                    residual -= entry.second * values[entry.first];
                  constraint_error = std::max(constraint_error, std::abs(residual));
                }
          }

        // A known coarse FE field must be prolonged without losing the rotation.
        // Use zero radial-boundary values, as required by the velocity correction.
        for (unsigned int level = 1; level < levels; ++level)
          {
            LevelVector coarse, prolonged, expected;
            coarse.reinit(exact[level-1]);
            coarse = exact[level-1];
            for (const auto index : dofs.locally_owned_mg_dofs(level-1))
              if (mg.is_boundary_index(level-1, index))
                coarse[index] = 0.;
            coarse.update_ghost_values();
            LevelVector reference = coarse;
            reference.update_ghost_values();
            const LevelVector &coarse_field = reference;
            coarse.zero_out_ghost_values();
            constraints[level-1].set_zero(coarse);
            prolonged.reinit(exact[level]);
            expected.reinit(exact[level]);
            transfer.prolongate(level, prolonged, coarse);
            prolonged.update_ghost_values();
            constraints[level].distribute(prolonged);

            std::vector<types::global_dof_index> parent_indices(fe.n_dofs_per_cell()), child_indices(fe.n_dofs_per_cell());
            for (const auto &cell : dofs.mg_cell_iterators_on_level(level))
              if (!cell->is_artificial_on_level())
                {
                  const auto parent = cell->parent();
                  parent->get_mg_dof_indices(parent_indices);
                  cell->get_mg_dof_indices(child_indices);
                  unsigned int child = 0;
                  while (parent->child(child) != cell)
                    ++child;
                  const auto &interpolation = fe.get_prolongation_matrix(child);
                  for (unsigned int i = 0; i < child_indices.size(); ++i)
                    if (dofs.locally_owned_mg_dofs(level).is_element(child_indices[i]))
                      {
                        double value = 0.;
                        for (unsigned int j = 0; j < parent_indices.size(); ++j)
                          value += interpolation(i,j) * coarse_field[parent_indices[j]];
                        expected[child_indices[i]] = value;
                      }
                }
            expected.compress(VectorOperation::insert);
            prolonged -= expected;
            transfer_error = std::max(transfer_error, prolonged.linfty_norm());
          }

        constraint_error = Utilities::MPI::max(constraint_error, this->get_mpi_communicator());
        transfer_error = Utilities::MPI::max(transfer_error, this->get_mpi_communicator());
        write_velocity_samples();
        this->get_pcout() << "   Velocity level constraint error: " << constraint_error << '\n'
                          << "   Velocity multigrid transfer error: " << transfer_error << std::endl;
        AssertThrow(constraint_error < 1.e-11, ExcMessage("Multigrid constraints do not preserve rigid rotation."));
        AssertThrow(transfer_error < 1.e-11, ExcMessage("Multigrid transfer does not preserve the rotated coarse FE field."));
        return {"Periodic velocity multigrid check:", "passed"};
      }

    private:
      void write_velocity_samples() const
      {
        FEFaceValues<dim> face_values(this->get_mapping(), this->get_fe(), QMidpoint<dim-1>(), update_values | update_quadrature_points);
        std::vector<Tensor<1,dim>> velocities(1);
        std::vector<std::array<double,6>> samples;
        const auto &material = Plugins::get_plugin_as_type<const AnnulusBenchmark::AnnulusMaterial<dim>>(this->get_material_model());
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned())
            for (const auto face : cell->face_indices())
              if (cell->has_periodic_neighbor(face))
                {
                  face_values.reinit(cell, face);
                  const auto point = face_values.quadrature_point(0);
                  if (point.norm() < 1.2 || point.norm() > 1.8)
                    continue;
                  face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), velocities);
                  const auto analytical = AnnulusBenchmark::AnalyticSolutions::Annulus_velocity(point, material.get_k(), this->get_time(), false);
                  samples.push_back({{point[0], point[1], velocities[0][0], velocities[0][1], analytical[0], analytical[1]}});
                }
        const auto collected = Utilities::MPI::gather(this->get_mpi_communicator(), samples);
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
            samples.clear();
            for (const auto &rank_samples : collected)
              samples.insert(samples.end(), rank_samples.begin(), rank_samples.end());
            AssertThrow(!samples.empty(), ExcMessage("No periodic velocity samples were written."));
            std::sort(samples.begin(), samples.end());
            std::ofstream output(this->get_output_directory() + "periodic_velocity.txt");
            output << "# x y ux uy analytical_ux analytical_uy\n" << std::scientific << std::setprecision(12);
            for (const auto &sample : samples)
              {
                for (const auto value : sample)
                  output << value << ' ';
                output << '\n';
              }
          }
      }
  };

  ASPECT_REGISTER_POSTPROCESSOR(PeriodicMultigridCheck,
                                "periodic multigrid check",
                                "Check the solver's velocity level constraints and transfer against a rotating FE field.")
}
