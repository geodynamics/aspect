/*
  Copyright (C) 2018 - 2024 by the authors of the ASPECT code.

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


#include <aspect/simulator/solver/stokes_matrix_free.h>
#include <aspect/simulator/solver/matrix_free_operators.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/melt.h>
#include <aspect/newton.h>

#include <deal.II/base/signaling_nan.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/read_write_vector.templates.h>
#include <deal.II/lac/solver_idr.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>

#include <deal.II/grid/manifold.h>

#include <deal.II/matrix_free/tools.h>

namespace aspect
{
  namespace internal
  {
    /**
     * Matrix-free operators must use deal.II defined vectors, while the rest of the ASPECT
     * software is based on Trilinos vectors. Here we define functions which copy between the
     * vector types.
     */
    namespace ChangeVectorTypes
    {
      void copy(aspect::LinearAlgebra::Vector &out,
                const dealii::LinearAlgebra::distributed::Vector<double> &in)
      {
        dealii::LinearAlgebra::ReadWriteVector<double> rwv(out.locally_owned_elements());
        rwv.import_elements(in, VectorOperation::insert);
        out.import_elements(rwv,VectorOperation::insert);
      }

      void copy(dealii::LinearAlgebra::distributed::Vector<double> &out,
                const aspect::LinearAlgebra::Vector &in)
      {
        dealii::LinearAlgebra::ReadWriteVector<double> rwv;
        rwv.reinit(in);
        out.import_elements(rwv, VectorOperation::insert);
      }

      void copy(aspect::LinearAlgebra::BlockVector &out,
                const dealii::LinearAlgebra::distributed::BlockVector<double> &in)
      {
        const unsigned int n_blocks = in.n_blocks();
        for (unsigned int b=0; b<n_blocks; ++b)
          copy(out.block(b),in.block(b));
      }

      void copy(dealii::LinearAlgebra::distributed::BlockVector<double> &out,
                const aspect::LinearAlgebra::BlockVector &in)
      {
        const unsigned int n_blocks = in.n_blocks();
        for (unsigned int b=0; b<n_blocks; ++b)
          copy(out.block(b),in.block(b));
      }
    }

  }


  namespace MatrixFreeStokesOperators
  {
    template <int dim, typename number>
    inline std::size_t
    OperatorCellData<dim,number>::memory_consumption() const
    {
      return viscosity.memory_consumption()
             + newton_factor_wrt_pressure_table.memory_consumption()
             + strain_rate_table.memory_consumption()
             + newton_factor_wrt_strain_rate_table.memory_consumption()
             + dilation_derivative_wrt_pressure_table.memory_consumption()
             + dilation_derivative_wrt_strain_rate_table.memory_consumption();
    }



    template <int dim, typename number>
    void
    OperatorCellData<dim,number>::clear()
    {
      enable_newton_derivatives = false;
      viscosity.clear();
      newton_factor_wrt_pressure_table.clear();
      strain_rate_table.clear();
      newton_factor_wrt_strain_rate_table.clear();
      dilation_derivative_wrt_pressure_table.clear();
      dilation_derivative_wrt_strain_rate_table.clear();
    }
  }



  template <int dim, int degree_v, typename number>
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::StokesOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>()
  {}



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::clear ()
  {
    this->cell_data = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::BlockVector<number>>::clear();
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::
  set_cell_data (const OperatorCellData<dim,number> &data)
  {
    this->cell_data = &data;
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>
  ::compute_diagonal ()
  {
    // There is currently no need in the code for the diagonal of the entire Stokes
    // block. If needed, one could easily construct based on the diagonal of the A
    // block and append zeros to the end for the number of pressure DoFs.
    Assert(false, ExcNotImplemented());
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>
  ::local_apply (const dealii::MatrixFree<dim, number>                         &data,
                 dealii::LinearAlgebra::distributed::BlockVector<number>       &dst,
                 const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                 const std::pair<unsigned int, unsigned int>                   &cell_range) const
  {
    FEEvaluation<dim,degree_v,degree_v+1,dim,number> u_eval(data, 0);
    FEEvaluation<dim,degree_v-1,degree_v+1,1,number> p_eval(data, /*dofh*/1);

    const bool use_viscosity_at_quadrature_points
      = (cell_data->viscosity.size(1) == u_eval.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        VectorizedArray<number> viscosity_x_2 = 2. * cell_data->viscosity(cell,0);

        u_eval.reinit(cell);
        u_eval.gather_evaluate(src.block(0), EvaluationFlags::gradients);

        p_eval.reinit(cell);
        p_eval.gather_evaluate(src.block(1), EvaluationFlags::values);

        // Derivative terms related to the Newton solver
        VectorizedArray<number> deta_deps_times_sym_grad_u(0.);
        VectorizedArray<number> eps_times_sym_grad_u(0.);
        VectorizedArray<number> deta_dp_times_p(0.);
        if (cell_data->enable_newton_derivatives)
          {
            SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u;
            VectorizedArray<number> val_p;
            for (const unsigned int q : u_eval.quadrature_point_indices())
              {
                sym_grad_u = u_eval.get_symmetric_gradient(q);
                val_p      = p_eval.get_value(q);
                deta_deps_times_sym_grad_u += cell_data->newton_factor_wrt_strain_rate_table(cell,q)
                                              * sym_grad_u;
                deta_dp_times_p += cell_data->newton_factor_wrt_pressure_table(cell,q) * val_p;
                if (cell_data->symmetrize_newton_system)
                  eps_times_sym_grad_u += cell_data->strain_rate_table(cell,q) * sym_grad_u;
              }
          }

        for (const unsigned int q : u_eval.quadrature_point_indices())
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              viscosity_x_2 = 2. * cell_data->viscosity(cell,q);

            const SymmetricTensor<2,dim,VectorizedArray<number>>
            sym_grad_u = u_eval.get_symmetric_gradient(q);
            const VectorizedArray<number> div_u = trace(sym_grad_u);
            const VectorizedArray<number> val_p = p_eval.get_value(q);

            // Terms to be tested by phi_p:
            VectorizedArray<number> pressure_terms =
              -cell_data->pressure_scaling * div_u;

            if (cell_data->enable_prescribed_dilation)
              pressure_terms -= cell_data->pressure_scaling *
                                cell_data->pressure_scaling *
                                cell_data->dilation_lhs_term_table(cell,q) *
                                val_p;

            // Terms to be tested by the symmetric gradients of phi_u:
            SymmetricTensor<2,dim,VectorizedArray<number>>
            velocity_terms = viscosity_x_2 * sym_grad_u;

            for (unsigned int d=0; d<dim; ++d)
              velocity_terms[d][d] -= cell_data->pressure_scaling * val_p;

            if (cell_data->is_compressible ||
                cell_data->enable_prescribed_dilation)
              for (unsigned int d=0; d<dim; ++d)
                velocity_terms[d][d] -= viscosity_x_2 / 3. * div_u;

            // Add the Newton derivatives if required.
            if (cell_data->enable_newton_derivatives)
              {
                velocity_terms +=
                  ( cell_data->symmetrize_newton_system ?
                    ( cell_data->strain_rate_table(cell,q) * deta_deps_times_sym_grad_u +
                      cell_data->newton_factor_wrt_strain_rate_table(cell,q) * eps_times_sym_grad_u ) :
                    2. * cell_data->strain_rate_table(cell,q) * deta_deps_times_sym_grad_u )
                  +
                  2. * cell_data->strain_rate_table(cell,q) * deta_dp_times_p;

                if (cell_data->enable_prescribed_dilation)
                  {
                    pressure_terms += ( ( cell_data->dilation_derivative_wrt_strain_rate_table(cell,q)
                                          * sym_grad_u )
                                        +
                                        ( cell_data->dilation_derivative_wrt_pressure_table(cell,q)
                                          * cell_data->pressure_scaling * val_p )
                                      )
                                      * cell_data->pressure_scaling;
                  }
              }

            u_eval.submit_symmetric_gradient(velocity_terms, q);
            p_eval.submit_value(pressure_terms, q);
          }

        u_eval.integrate_scatter(EvaluationFlags::gradients, dst.block(0));
        p_eval.integrate_scatter(EvaluationFlags::values, dst.block(1));
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim, degree_v, number>
  ::local_apply_face(const dealii::MatrixFree<dim, number> &,
                     dealii::LinearAlgebra::distributed::BlockVector<number> &,
                     const dealii::LinearAlgebra::distributed::BlockVector<number> &,
                     const std::pair<unsigned int, unsigned int> &) const
  {
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim, degree_v, number>
  ::local_apply_boundary_face(const dealii::MatrixFree<dim, number> &data,
                              dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                              const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                              const std::pair<unsigned int, unsigned int> &face_range) const
  {
    // Assemble the fictive stabilization stress (phi_u[i].g)*(phi_u[j].n)
    // g=pressure_perturbation * g_hat is stored in free_surface_stabilization_term_table
    //  n is the normal vector
    FEFaceEvaluation<dim, degree_v, degree_v + 1, dim, number> velocity(data);
    const unsigned int n_faces_interior = data.n_inner_face_batches();

    for (unsigned int face = face_range.first; face < face_range.second; ++face)
      {
        const auto boundary_id = data.get_boundary_id(face);
        if (cell_data->free_surface_boundary_indicators.find(boundary_id)
            == cell_data->free_surface_boundary_indicators.end())
          continue;

        velocity.reinit(face);
        velocity.gather_evaluate (src.block(0), EvaluationFlags::values);

        for (const unsigned int q : velocity.quadrature_point_indices())
          {
            const Tensor<1, dim, VectorizedArray<number>> phi_u_i = velocity.get_value(q);
#if DEAL_II_VERSION_GTE(9,7,0)
            const auto &normal_vector = velocity.normal_vector(q);
#else
            const auto &normal_vector = velocity.get_normal_vector(q);
#endif
            const auto stabilization_tensor = cell_data->free_surface_stabilization_term_table(face - n_faces_interior, q);
            const auto value_submit = -(stabilization_tensor * phi_u_i) * normal_vector;

            velocity.submit_value(value_submit, q);
          }
        velocity.integrate_scatter(EvaluationFlags::values, dst.block(0));
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>
  ::apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
               const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const
  {
    if (cell_data->apply_stabilization_free_surface_faces)
      MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>::
      data->loop(&StokesOperator::local_apply,
                 &StokesOperator::local_apply_face,
                 &StokesOperator::local_apply_boundary_face,
                 this,
                 dst,
                 src,
                 false, /*zero_dst_vector*/
                 MatrixFree<dim, number>::DataAccessOnFaces::values,
                 MatrixFree<dim, number>::DataAccessOnFaces::values);

    else
      MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>::
      data->cell_loop(&StokesOperator::local_apply, this, dst, src);
  }

  template <int dim, int degree_v, typename number>
  MatrixFreeStokesOperators::BTBlockOperator<dim,degree_v,number>::BTBlockOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>()
  {}



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::BTBlockOperator<dim,degree_v,number>::clear ()
  {
    this->cell_data = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::BlockVector<number>>::clear();
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::BTBlockOperator<dim,degree_v,number>::
  set_cell_data (const OperatorCellData<dim,number> &data)
  {
    this->cell_data = &data;
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::BTBlockOperator<dim,degree_v,number>
  ::compute_diagonal ()
  {
    // There no need in the code for this diagonal.
    Assert(false, ExcNotImplemented());
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::BTBlockOperator<dim,degree_v,number>
  ::local_apply (const dealii::MatrixFree<dim, number>                         &data,
                 dealii::LinearAlgebra::distributed::BlockVector<number>       &dst,
                 const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                 const std::pair<unsigned int, unsigned int>                   &cell_range) const
  {
    FEEvaluation<dim,degree_v,degree_v+1,dim,number> u_eval(data, 0);
    FEEvaluation<dim,degree_v-1,degree_v+1,1,number> p_eval(data, /*dofh*/1);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        u_eval.reinit(cell);

        p_eval.reinit(cell);
        p_eval.gather_evaluate(src.block(1), EvaluationFlags::values);

        for (const unsigned int q : u_eval.quadrature_point_indices())
          {
            const VectorizedArray<number> val_p = p_eval.get_value(q);


            SymmetricTensor<2,dim,VectorizedArray<number>>
            velocity_terms;

            for (unsigned int d=0; d<dim; ++d)
              velocity_terms[d][d] -= cell_data->pressure_scaling * val_p;
            u_eval.submit_symmetric_gradient(velocity_terms, q);
          }

        u_eval.integrate_scatter(EvaluationFlags::gradients, dst.block(0));
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::BTBlockOperator<dim, degree_v, number>
  ::local_apply_face(const dealii::MatrixFree<dim, number> &,
                     dealii::LinearAlgebra::distributed::BlockVector<number> &,
                     const dealii::LinearAlgebra::distributed::BlockVector<number> &,
                     const std::pair<unsigned int, unsigned int> &) const
  {
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::BTBlockOperator<dim,degree_v,number>
  ::apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
               const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const
  {
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number>>::
    data->cell_loop(&BTBlockOperator::local_apply, this, dst, src);
  }

  /**
   * Mass matrix operator on pressure
   */
  template <int dim, int degree_p, typename number>
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::MassMatrixOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>()
  {}



  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::clear ()
  {
    this->cell_data = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number>>::clear();
  }



  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::reinit(const Mapping<dim>              &mapping,
                                                                             const DoFHandler<dim>           &dof_handler_v,
                                                                             const DoFHandler<dim>           &dof_handler_p,
                                                                             const AffineConstraints<number> &constraints_v,
                                                                             const AffineConstraints<number> &constraints_p,
                                                                             std::shared_ptr<MatrixFree<dim,double>> mf_storage,
                                                                             const unsigned int level)
  {
    typename MatrixFree<dim, number>::AdditionalData data;
    data.mapping_update_flags =
      update_quadrature_points /*| update_gradients*/ | update_values;
    data.mg_level = level;

    data.tasks_parallel_scheme =
      MatrixFree<dim,double>::AdditionalData::none;
    AffineConstraints<number> dummy;

    mf_storage->reinit(mapping,
    std::vector< const DoFHandler< dim > *> {&dof_handler_v, &dof_handler_p},
    std::vector< const AffineConstraints< number > *> {&constraints_v, &constraints_p} ,
    QGauss<1>(degree_p+2), data);

    this->initialize(mf_storage, std::vector< unsigned int > {1}, std::vector< unsigned int > {1});
  }



  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::
  set_cell_data (const OperatorCellData<dim,number> &data)
  {
    this->cell_data = &data;
  }



  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                 dealii::LinearAlgebra::distributed::Vector<number>       &dst,
                 const dealii::LinearAlgebra::distributed::Vector<number> &src,
                 const std::pair<unsigned int, unsigned int>           &cell_range) const
  {
    FEEvaluation<dim,degree_p,degree_p+2,1,number> pressure (data, /*dofh*/1);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        pressure.reinit (cell);
        pressure.gather_evaluate (src, EvaluationFlags::values);
        this->inner_cell_operation(pressure);
        pressure.integrate_scatter (EvaluationFlags::values, dst);
      }

  }


  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::inner_cell_operation(FEEvaluation<dim,
                         degree_p,
                         degree_p+2,
                         1,
                         number> &pressure) const
  {
    const bool use_viscosity_at_quadrature_points
      = (cell_data->viscosity.size(1) == pressure.n_q_points);

    const unsigned int cell = pressure.get_current_cell_index();
    const unsigned int n_components_filled = this->get_matrix_free()->n_active_entries_per_cell_batch(cell);

    VectorizedArray<number> prefactor;

    // The /= operator for VectorizedArray results in a floating point operation
    // (divide by 0) since the (*viscosity)(cell) array is not completely filled.
    // Therefore, we need to divide each entry manually.
    if (!use_viscosity_at_quadrature_points)
      {
        for (unsigned int c=0; c<n_components_filled; ++c)
          prefactor[c] = cell_data->pressure_scaling*cell_data->pressure_scaling / cell_data->viscosity(cell, 0)[c];
      }

    for (const unsigned int q : pressure.quadrature_point_indices())
      {
        // Only update the viscosity if a Q1 projection is used.
        if (use_viscosity_at_quadrature_points)
          {
            for (unsigned int c=0; c<n_components_filled; ++c)
              prefactor[c] = cell_data->pressure_scaling*cell_data->pressure_scaling / cell_data->viscosity(cell, q)[c];
          }

        pressure.submit_value(prefactor*pressure.get_value(q), q);
      }
  }



  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
               const dealii::LinearAlgebra::distributed::Vector<number> &src) const
  {
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number>>::
    data->cell_loop(&MassMatrixOperator::local_apply, this, dst, src);
  }



  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::compute_diagonal ()
  {
    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number>>>();
    this->diagonal_entries =
      std::make_shared<DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number>>>();

    dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    dealii::LinearAlgebra::distributed::Vector<number> &diagonal =
      this->diagonal_entries->get_vector();

    this->data->initialize_dof_vector(inverse_diagonal, /*dofh*/1);
    this->data->initialize_dof_vector(diagonal, /*dofh*/1);

    MatrixFreeTools::compute_diagonal<dim,degree_p,degree_p+2,1,number,VectorizedArray<number>,dealii::LinearAlgebra::distributed::Vector<number>>(
      *(this->get_matrix_free()),
      diagonal,
      [&](FEEvaluation<dim,
          degree_p,
          degree_p+2,
          1,
          number> &pressure)
    {
      pressure.evaluate(EvaluationFlags::values);
      this->inner_cell_operation(pressure);
      pressure.integrate(EvaluationFlags::values);
    },
    1 /* dofhandler */);

    this->set_constrained_entries_to_one(diagonal);
    inverse_diagonal = diagonal;

    // Finally loop over all of the computed diagonal elements and invert them.
    // The following loop relies on the fact that inverse_diagonal.begin()/end()
    // iterates only over the *locally owned* elements of the vector in which
    // we store inverse_diagonal.
    for (auto &local_element : inverse_diagonal)
      {
        Assert(local_element > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero or negative."));
        local_element = 1./local_element;
      }
  }



  /**
   * Velocity block operator
   */
  template <int dim, int degree_v, typename number>
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::ABlockOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>()
  {}



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::clear ()
  {
    this->cell_data = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number>>::clear();
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::reinit(const Mapping<dim>              &mapping,
                                                                         const DoFHandler<dim>           &dof_handler_v,
                                                                         const DoFHandler<dim>           &dof_handler_p,
                                                                         const AffineConstraints<number> &constraints_v,
                                                                         const AffineConstraints<number> &constraints_p,
                                                                         std::shared_ptr<MatrixFree<dim,double>> mf_storage,
                                                                         const unsigned int level)
  {
    typename MatrixFree<dim, number>::AdditionalData data;
    data.mapping_update_flags = update_quadrature_points | update_gradients | update_values;
    data.mg_level = level;

    typename MatrixFree<dim,double>::AdditionalData additional_data;
    additional_data.tasks_parallel_scheme =
      MatrixFree<dim,double>::AdditionalData::none;
    additional_data.mapping_update_flags = (update_values | update_gradients |
                                            update_JxW_values | update_quadrature_points);

    AffineConstraints<number> dummy;
    mf_storage->reinit(mapping,
    std::vector< const DoFHandler< dim > *> {&dof_handler_v, &dof_handler_p},
    std::vector< const AffineConstraints< number > *> {&constraints_v, &constraints_p} ,
    QGauss<1>(degree_v+1), additional_data);

    this->initialize(mf_storage, std::vector< unsigned int > {0}, std::vector< unsigned int > {0});
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::
  set_cell_data (const OperatorCellData<dim,number> &data)
  {
    this->cell_data = &data;
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::inner_cell_operation(FEEvaluation<dim,
                         degree_v,
                         degree_v+1,
                         dim,
                         number> &velocity) const
  {
    const bool use_viscosity_at_quadrature_points
      = (cell_data->viscosity.size(1) == velocity.n_q_points);

    const unsigned int cell = velocity.get_current_cell_index();
    VectorizedArray<number> viscosity_x_2 = 2.0*cell_data->viscosity(cell, 0);

    for (const unsigned int q : velocity.quadrature_point_indices())
      {
        // Only update the viscosity if a Q1 projection is used.
        if (use_viscosity_at_quadrature_points)
          viscosity_x_2 = 2.0*cell_data->viscosity(cell, q);

        SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
          velocity.get_symmetric_gradient (q);
        sym_grad_u *= viscosity_x_2;

        if (cell_data->is_compressible ||
            cell_data->enable_prescribed_dilation)
          {
            const VectorizedArray<number> div = trace(sym_grad_u);
            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= 1.0/3.0*div;
          }

        velocity.submit_symmetric_gradient(sym_grad_u, q);
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::cell_operation(FEEvaluation<dim,
                   degree_v,
                   degree_v+1,
                   dim,
                   number> &velocity) const
  {
    velocity.evaluate (EvaluationFlags::gradients);
    this->inner_cell_operation(velocity);
    velocity.integrate(EvaluationFlags::gradients);
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                 dealii::LinearAlgebra::distributed::Vector<number>       &dst,
                 const dealii::LinearAlgebra::distributed::Vector<number> &src,
                 const std::pair<unsigned int, unsigned int>           &cell_range) const
  {
    FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (data,0);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        velocity.reinit (cell);

        // Instead of calling
        //   velocity.read_dof_values (src);
        //   velocity.evaluate (EvaluationFlags::gradients);
        // (the latter by calling cell_operation()), we use the more efficient
        // combined gather_evaluate() and use inner_cell_operation().
        velocity.gather_evaluate (src, EvaluationFlags::gradients);
        this->inner_cell_operation(velocity);
        velocity.integrate_scatter (EvaluationFlags::gradients, dst);
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
               const dealii::LinearAlgebra::distributed::Vector<number> &src) const
  {
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number>>::
    data->cell_loop(&ABlockOperator::local_apply, this, dst, src);
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::compute_diagonal ()
  {
    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number>>>();
    dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);

    MatrixFreeTools::compute_diagonal(
      *(this->get_matrix_free()),
      inverse_diagonal,
      &MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::cell_operation,
      this);

    this->set_constrained_entries_to_one(inverse_diagonal);

    // Finally loop over all of the computed diagonal elements and invert them.
    // The following loop relies on the fact that inverse_diagonal.begin()/end()
    // iterates only over the *locally owned* elements of the vector in which
    // we store inverse_diagonal.
    for (auto &local_element : inverse_diagonal)
      {
        Assert(local_element > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero or negative."));
        local_element = 1./local_element;
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::set_diagonal (const dealii::LinearAlgebra::distributed::Vector<number> &diag)
  {
    this->inverse_diagonal_entries =
      std::make_shared<DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number>>>();
    dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);

    inverse_diagonal = diag;

    this->set_constrained_entries_to_one(inverse_diagonal);

    // Finally loop over all of the computed diagonal elements and invert them.
    // The following loop relies on the fact that inverse_diagonal.begin()/end()
    // iterates only over the *locally owned* elements of the vector in which
    // we store inverse_diagonal.
    for (auto &local_element : inverse_diagonal)
      {
        Assert(local_element > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero or negative."));
        local_element = 1./local_element;
      }
  }

}

// explicit instantiations
namespace aspect
{
#define INSTANTIATE(dim) \
  template class MatrixFreeStokesOperators::ABlockOperator<dim,2,GMGNumberType>; \
  template class MatrixFreeStokesOperators::ABlockOperator<dim,3,GMGNumberType>; \
  template class MatrixFreeStokesOperators::StokesOperator<dim,2,GMGNumberType>; \
  template class MatrixFreeStokesOperators::StokesOperator<dim,3,GMGNumberType>; \
  template class MatrixFreeStokesOperators::BTBlockOperator<dim,2,GMGNumberType>; \
  template class MatrixFreeStokesOperators::BTBlockOperator<dim,3,GMGNumberType>; \
  template class MatrixFreeStokesOperators::MassMatrixOperator<dim,1,GMGNumberType>; \
  template class MatrixFreeStokesOperators::MassMatrixOperator<dim,2,GMGNumberType>; \
  template struct MatrixFreeStokesOperators::OperatorCellData<dim, GMGNumberType>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE

}
