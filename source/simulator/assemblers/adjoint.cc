/*
  Copyright (C) 2016 - 2020 by the authors of the ASPECT code.

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

#include <aspect/simulator.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <aspect/simulator/assemblers/adjoint.h>
#include <aspect/postprocess/dynamic_topography.h>

namespace aspect
{
  namespace Assemblers
  {

    template <int dim>
    void
    StokesAdjointRHS<dim>::
    execute (internal::Assembly::Scratch::ScratchBase<dim>   &scratch_base,
             internal::Assembly::CopyData::CopyDataBase<dim> &data_base) const
    {
      internal::Assembly::Scratch::StokesSystem<dim> &scratch = dynamic_cast<internal::Assembly::Scratch::StokesSystem<dim>& > (scratch_base);
      internal::Assembly::CopyData::StokesSystem<dim> &data = dynamic_cast<internal::Assembly::CopyData::StokesSystem<dim>& > (data_base);

      const Introspection<dim> &introspection = this->introspection();
      const FiniteElement<dim> &fe = scratch.finite_element_values.get_fe();

      const unsigned int stokes_dofs_per_cell = data.local_dof_indices.size();
      const unsigned int n_face_q_points      = scratch.face_finite_element_values.n_quadrature_points;

      const double pressure_scaling = this->get_pressure_scaling();

      // TODO: make this consistent with input eventually
      const double density_above = 0.;

      // Get a pointer to the dynamic topography postprocessor.
      const Postprocess::DynamicTopography<dim> &dynamic_topography =
        this->get_postprocess_manager().template get_matching_postprocessor<Postprocess::DynamicTopography<dim> >();

      // Get the already-computed dynamic topography solution.
      const LinearAlgebra::BlockVector &topography_vector = dynamic_topography.topography_vector();
      std::vector<double> topo_values( n_face_q_points );

      // check that the cell is at the top and that the cell is at the top
      if (scratch.cell->face(scratch.face_number)->at_boundary()
          &&
          this->get_geometry_model().depth (scratch.cell->face(scratch.face_number)->center()) < scratch.cell->face(scratch.face_number)->minimum_vertex_distance()/3)
        {
          // get values at the surface of the cell
          scratch.face_finite_element_values.reinit (scratch.cell, scratch.face_number);
          scratch.face_finite_element_values[introspection.extractors.temperature].get_function_values(topography_vector, topo_values);

          // ----------  Assemble RHS  ---------------

          for (unsigned int q=0; q<n_face_q_points; ++q)
            {
              for (unsigned int i=0, i_stokes=0; i_stokes<stokes_dofs_per_cell; /*increment at end of loop*/)
                {
                  if (introspection.is_stokes_component(fe.system_to_component_index(i).first))
                    {
                      scratch.phi_p[i_stokes] = scratch.face_finite_element_values[introspection.extractors.pressure].value (i, q);
                      scratch.grads_phi_u[i_stokes] = scratch.face_finite_element_values[introspection.extractors.velocities].symmetric_gradient(i,q);
                      ++i_stokes;
                    }
                  ++i;
                }


              const Tensor<1,dim> n_hat = scratch.face_finite_element_values.normal_vector(q);
              const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector (scratch.face_finite_element_values.quadrature_point(q));
              const double density = scratch.material_model_outputs.densities[q];
              const double eta = scratch.material_model_outputs.viscosities[q];
              const double JxW = scratch.face_finite_element_values.JxW(q);

              for (unsigned int i=0; i<stokes_dofs_per_cell; ++i)
                {
                  data.local_rhs(i) += topo_values[q] * (2.0*eta *(n_hat * (scratch.grads_phi_u[i] * n_hat))
                                                           - pressure_scaling *scratch.phi_p[i]) / ((density-density_above)* gravity.norm())
                                       * JxW;
                }
            }
        }
    }

  }

} // namespace aspect


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Assemblers
  {

#define INSTANTIATE(dim) \
  template class StokesAdjointRHS<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE

  }
}
