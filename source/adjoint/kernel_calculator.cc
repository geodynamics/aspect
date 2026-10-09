/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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

#include <aspect/adjoint/kernel_calculator.h>
#include <aspect/simulator.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <limits>

namespace aspect
{
  namespace Adjoint
  {
    template <int dim>
    void
    KernelRepository<dim>::add_contribution(const KernelContributionKey &key,
                                            const Vector<double> &values)
    {
      AssertThrow(contribution_map.find(key) == contribution_map.end(),
                  ExcMessage("Duplicate adjoint kernel contribution for objective <"
                             + key.objective_name + ">, physics term <"
                             + key.physics_term_name + ">, and physical property <"
                             + property_name(key.property) + ">."));
      contribution_map[key] = values;
    }



    template <int dim>
    void
    KernelRepository<dim>::set_cell_volumes(const Vector<double> &values)
    {
      cell_volume_values = values;
    }



    template <int dim>
    const Vector<double> &
    KernelRepository<dim>::cell_volumes() const
    {
      return cell_volume_values;
    }



    template <int dim>
    bool
    KernelRepository<dim>::empty() const
    {
      return contribution_map.empty();
    }



    template <int dim>
    unsigned int
    KernelRepository<dim>::n_contributions() const
    {
      return contribution_map.size();
    }



    template <int dim>
    const std::map<KernelContributionKey, Vector<double>> &
    KernelRepository<dim>::contributions() const
    {
      return contribution_map;
    }



    template <int dim>
    KernelCalculator<dim>::KernelCalculator(const Simulator<dim> &simulator)
    {
      this->initialize_simulator(simulator);
    }



    template <int dim>
    KernelRepository<dim>
    KernelCalculator<dim>::calculate(
      const std::vector<std::unique_ptr<ObjectiveResult<dim>>> &objective_results,
      const std::vector<std::unique_ptr<AdjointState<dim>>> &adjoint_states) const
    {
      AssertThrow(objective_results.size() == adjoint_states.size(),
                  ExcMessage("The number of objective results must match the number of adjoint states."));

      KernelRepository<dim> kernels;

      const unsigned int n_active_cells = this->get_triangulation().n_active_cells();
      const unsigned int quadrature_degree =
        this->get_fe().base_element(this->introspection().base_elements.velocities).degree + 1;
      const QGauss<dim> quadrature_formula(quadrature_degree);

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature_formula,
                              update_values |
                              update_gradients |
                              update_quadrature_points |
                              update_JxW_values);

      for (unsigned int objective_index = 0; objective_index < objective_results.size(); ++objective_index)
        {
          const std::string &objective_name = objective_results[objective_index]->objective_name;
          const LinearAlgebra::BlockVector &adjoint_solution = adjoint_states[objective_index]->solution;

          Vector<double> cell_volumes(n_active_cells);
          Vector<double> density_volume(n_active_cells);
          Vector<double> viscosity_volume(n_active_cells);

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                fe_values.reinit(cell);

                MaterialModel::MaterialModelInputs<dim> adjoint_in(fe_values,
                                                                   cell,
                                                                   this->introspection(),
                                                                   adjoint_solution);

                std::vector<double> local_forward_values(this->get_fe().dofs_per_cell);
                std::vector<double> local_adjoint_values(this->get_fe().dofs_per_cell);
                cell->get_dof_values(this->get_solution(),
                                     local_forward_values.begin(),
                                     local_forward_values.end());
                cell->get_dof_values(adjoint_solution,
                                     local_adjoint_values.begin(),
                                     local_adjoint_values.end());

                const unsigned int cell_index = cell->active_cell_index();
                for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
                  {
                    const Tensor<1,dim> gravity =
                      this->get_gravity_model().gravity_vector(fe_values.quadrature_point(q));
                    const double JxW = fe_values.JxW(q);
                    SymmetricTensor<2, dim> forward_operator_strain;
                    SymmetricTensor<2, dim> adjoint_operator_strain;

                    for (unsigned int i = 0; i < this->get_fe().dofs_per_cell; ++i)
                      if (this->introspection().is_stokes_component(this->get_fe().system_to_component_index(i).first))
                        {
                          const SymmetricTensor<2, dim> epsilon_phi_u =
                            fe_values[this->introspection().extractors.velocities].symmetric_gradient(i, q);

                          forward_operator_strain += local_forward_values[i] * epsilon_phi_u;
                          adjoint_operator_strain += local_adjoint_values[i] * epsilon_phi_u;
                        }

                    cell_volumes(cell_index) += JxW;
                    density_volume(cell_index) += -(gravity * adjoint_in.velocity[q]) * JxW;
                    // Physical-property viscosity controls are additive
                    // viscosity increments. The incompressible Stokes
                    // assembler uses the full symmetric gradient term
                    // 2 eta eps(u):eps(v), so the derivative is with respect
                    // to eta itself and does not include an extra eta factor
                    // or a deviatoric projection.
                    viscosity_volume(cell_index) +=
                      (2.0 * (forward_operator_strain * adjoint_operator_strain))
                      * JxW;

                  }

              }

          for (const auto &cell : this->get_dof_handler().active_cell_iterators())
            if (cell->is_locally_owned())
              {
                const unsigned int cell_index = cell->active_cell_index();
                AssertThrow(cell_volumes(cell_index) > std::numeric_limits<double>::min(),
                            ExcMessage("KernelCalculator encountered a locally owned cell with zero volume."));

                density_volume(cell_index) /= cell_volumes(cell_index);
                viscosity_volume(cell_index) /= cell_volumes(cell_index);
              }


          if (objective_index == 0)
            kernels.set_cell_volumes(cell_volumes);

          kernels.add_contribution({objective_name,
                                    "incompressible Stokes volume",
                                    PhysicalProperty::density
                                   },
                                   density_volume);
          kernels.add_contribution({objective_name,
                                    "incompressible Stokes volume",
                                    PhysicalProperty::viscosity
                                   },
                                   viscosity_volume);
        }

      return kernels;
    }

  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Adjoint::KernelRepository<dim>; \
  template class Adjoint::KernelCalculator<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
