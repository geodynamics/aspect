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


#include <aspect/stokes_matrix_free.h>
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
      void copy(TrilinosWrappers::MPI::Vector &out,
                const dealii::LinearAlgebra::distributed::Vector<double> &in)
      {
        dealii::LinearAlgebra::ReadWriteVector<double> rwv(out.locally_owned_elements());
        rwv.import_elements(in, VectorOperation::insert);
        out.import_elements(rwv,VectorOperation::insert);
      }

      void copy(dealii::LinearAlgebra::distributed::Vector<double> &out,
                const TrilinosWrappers::MPI::Vector &in)
      {
        dealii::LinearAlgebra::ReadWriteVector<double> rwv;
        rwv.reinit(in);
        out.import_elements(rwv, VectorOperation::insert);
      }

      void copy(TrilinosWrappers::MPI::BlockVector &out,
                const dealii::LinearAlgebra::distributed::BlockVector<double> &in)
      {
        const unsigned int n_blocks = in.n_blocks();
        for (unsigned int b=0; b<n_blocks; ++b)
          copy(out.block(b),in.block(b));
      }

      void copy(dealii::LinearAlgebra::distributed::BlockVector<double> &out,
                const TrilinosWrappers::MPI::BlockVector &in)
      {
        const unsigned int n_blocks = in.n_blocks();
        for (unsigned int b=0; b<n_blocks; ++b)
          copy(out.block(b),in.block(b));
      }
    }


    /**
     * Implement the block Schur preconditioner for the Stokes system.
     */
    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    class BlockSchurGMGPreconditioner : public Subscriptor
    {
      public:
        /**
         * @brief Constructor
         *
         * @param Stokes_matrix The entire Stokes matrix
         * @param A_block The A block of the Stokes matrix
         * @param Schur_complement_block The matrix which describes the Schur complement approximation
         * @param A_block_preconditioner Preconditioner object for the matrix A.
         * @param Schur_complement_preconditioner Preconditioner object for the Schur complement.
         * @param do_solve_A A flag indicating whether we should actually solve with
         *     the matrix $A_block$, or only apply one preconditioner step with it.
         * @param do_solve_Schur_complement A flag indicating whether we should actually solve with
         *     the matrix $Schur_complement_block$, or only apply one preconditioner step with it.
         * @param A_block_is_symmetric A flag indicating whether the A block is symmetric.
         * @param A_block_tolerance The tolerance for the CG solver which computes
         *     the inverse of the A block.
         * @param Schur_complement_tolerance The tolerance for the CG solver which computes
         *     the inverse of the Schur complement block (Schur complement approximation matrix).
         */
        BlockSchurGMGPreconditioner (const StokesMatrixType                  &Stokes_matrix,
                                     const ABlockMatrixType                  &A_block,
                                     const SchurComplementMatrixType         &Schur_complement_block,
                                     const ABlockPreconditionerType          &A_block_preconditioner,
                                     const SchurComplementPreconditionerType &Schur_complement_preconditioner,
                                     const bool                               do_solve_A,
                                     const bool                               do_solve_Schur_complement,
                                     const bool                               A_block_is_symmetric,
                                     const double                             A_block_tolerance,
                                     const double                             Schur_complement_tolerance);

        /**
         * Matrix vector product with this preconditioner object.
         */
        void vmult (dealii::LinearAlgebra::distributed::BlockVector<double>       &dst,
                    const dealii::LinearAlgebra::distributed::BlockVector<double> &src) const;

        unsigned int n_iterations_A_block() const;
        unsigned int n_iterations_Schur_complement() const;


      private:
        /**
         * References to the various matrix object this preconditioner works on.
         */
        const StokesMatrixType                  &stokes_matrix;
        const ABlockMatrixType                  &A_block;
        const SchurComplementMatrixType         &Schur_complement_block;
        const ABlockPreconditionerType          &A_block_preconditioner;
        const SchurComplementPreconditionerType &Schur_complement_preconditioner;

        /**
         * Whether to actually invert the $\tilde M$ or $\tilde A$ of the preconditioner matrix
         * or to just apply a single preconditioner step with it.
         */
        const bool                                                      do_solve_A;
        const bool                                                      do_solve_Schur_complement;
        const bool                                                      A_block_is_symmetric;
        mutable unsigned int                                            n_iterations_A_;
        mutable unsigned int                                            n_iterations_Schur_complement_;
        const double                                                    A_block_tolerance;
        const double                                                    Schur_complement_tolerance;
        mutable dealii::LinearAlgebra::distributed::BlockVector<double> utmp;
    };

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                BlockSchurGMGPreconditioner (const StokesMatrixType                  &Stokes_matrix,
                                                             const ABlockMatrixType                  &A_block,
                                                             const SchurComplementMatrixType         &Schur_complement_block,
                                                             const ABlockPreconditionerType          &A_block_preconditioner,
                                                             const SchurComplementPreconditionerType &Schur_complement_preconditioner,
                                                             const bool                               do_solve_A,
                                                             const bool                               do_solve_Schur_complement,
                                                             const bool                               A_block_symmetric,
                                                             const double                             A_block_tolerance,
                                                             const double                             Schur_complement_tolerance)
                                  :
                                  stokes_matrix                   (Stokes_matrix),
                                  A_block                         (A_block),
                                  Schur_complement_block          (Schur_complement_block),
                                  A_block_preconditioner          (A_block_preconditioner),
                                  Schur_complement_preconditioner (Schur_complement_preconditioner),
                                  do_solve_A                      (do_solve_A),
                                  do_solve_Schur_complement       (do_solve_Schur_complement),
                                  A_block_is_symmetric            (A_block_symmetric),
                                  n_iterations_A_                 (0),
                                  n_iterations_Schur_complement_  (0),
                                  A_block_tolerance               (A_block_tolerance),
                                  Schur_complement_tolerance      (Schur_complement_tolerance)
    {}

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    unsigned int
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                n_iterations_A_block() const
    {
      return n_iterations_A_;
    }

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    unsigned int
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                n_iterations_Schur_complement() const
    {
      return n_iterations_Schur_complement_;
    }

    template <class StokesMatrixType, class ABlockMatrixType, class SchurComplementMatrixType,
              class ABlockPreconditionerType, class SchurComplementPreconditionerType>
    void
    BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType,
                                ABlockPreconditionerType, SchurComplementPreconditionerType>::
                                vmult (dealii::LinearAlgebra::distributed::BlockVector<double>       &dst,
                                       const dealii::LinearAlgebra::distributed::BlockVector<double>  &src) const
    {
      if (utmp.size()==0)
        utmp.reinit(src);

      // This needs to be done explicitly, as GMRES does not
      // initialize the data of the vector dst before calling
      // us. Otherwise we might use random data as our initial guess.
      dst = 0.0;

      // either solve with the Schur complement matrix (if do_solve_Schur_complement==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_Schur_complement)
        {
          // first solve with the bottom right block, which we have built
          // as a mass matrix with the inverse of the viscosity
          SolverControl solver_control(100, src.block(1).l2_norm() * Schur_complement_tolerance,true);

          SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> solver(solver_control);
          // Trilinos reports a breakdown
          // in case src=dst=0, even
          // though it should return
          // convergence without
          // iterating. We simply skip
          // solving in this case.
          if (src.block(1).l2_norm() > 1e-50)
            {
              try
                {
                  solver.solve(Schur_complement_block,
                               dst.block(1), src.block(1),
                               Schur_complement_preconditioner);
                  n_iterations_Schur_complement_ += solver_control.last_step();
                }
              // if the solver fails, report the error from processor 0 with some additional
              // information about its location, and throw a quiet exception on all other
              // processors
              catch (const std::exception &exc)
                {
                  Utilities::throw_linear_solver_failure_exception("iterative (bottom right) solver",
                                                                   "BlockSchurGMGPreconditioner::vmult",
                                                                   std::vector<SolverControl> {solver_control},
                                                                   exc,
                                                                   src.block(0).get_mpi_communicator());
                }
            }
        }
      else
        {
          Schur_complement_preconditioner.vmult(dst.block(1),src.block(1));
          n_iterations_Schur_complement_ += 1;
        }

      dst.block(1) *= -1.0;

      {
        // the u-block of dst only contains zeros
        stokes_matrix.vmult(utmp, dst); // B^T
        utmp.block(0) *= -1.0;
        utmp.block(0) += src.block(0);
      }

      // now either solve with the top left block (if do_solve_A==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          SolverControl solver_control(1000, utmp.block(0).l2_norm() * A_block_tolerance);
          PrimitiveVectorMemory<dealii::LinearAlgebra::distributed::Vector<double>> mem;

          try
            {
              if (A_block_is_symmetric)
                {
                  SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> solver(solver_control,mem);
                  solver.solve(A_block, dst.block(0), utmp.block(0),
                               A_block_preconditioner);
                }
              else
                {
                  // Use BiCGStab for non-symmetric matrices.
                  // BiCGStab can also solve indefinite systems if necessary.
                  // Do not compute the exact residual, as this
                  // is more expensive, and we only need an approximate solution.
                  SolverBicgstab<dealii::LinearAlgebra::distributed::Vector<double>>
                  solver(solver_control,
                         mem,
                         SolverBicgstab<dealii::LinearAlgebra::distributed::Vector<double>>::AdditionalData(/*exact_residual=*/ false));
                  solver.solve(A_block, dst.block(0), utmp.block(0),
                               A_block_preconditioner);
                }

              n_iterations_A_ += solver_control.last_step();
            }
          // if the solver fails, report the error from processor 0 with some additional
          // information about its location, and throw a quiet exception on all other
          // processors
          catch (const std::exception &exc)
            {
              Utilities::throw_linear_solver_failure_exception("iterative (top left) solver",
                                                               "BlockSchurGMGPreconditioner::vmult",
                                                               std::vector<SolverControl> {solver_control},
                                                               exc,
                                                               src.block(0).get_mpi_communicator());
            }
        }
      else
        {
          A_block_preconditioner.vmult (dst.block(0), utmp.block(0));
          n_iterations_A_ += 1;
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
             + newton_factor_wrt_strain_rate_table.memory_consumption();
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

        // Store the symmetric gradients of the velocity field and the
        // values of the pressure field
        AlignedVector<SymmetricTensor<2,dim,VectorizedArray<number>>> sym_grad_u;
        AlignedVector<VectorizedArray<number>> val_p;
        if (cell_data->enable_newton_derivatives)
          {
            sym_grad_u.resize(u_eval.n_q_points);
            val_p.resize(u_eval.n_q_points);
            for (const unsigned int r : u_eval.quadrature_point_indices())
              {
                sym_grad_u[r] = u_eval.get_symmetric_gradient(r);
                val_p[r]      = p_eval.get_value(r);
              }
          }

        for (const unsigned int q : u_eval.quadrature_point_indices())
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              viscosity_x_2 = 2. * cell_data->viscosity(cell,q);

            const SymmetricTensor<2,dim,VectorizedArray<number>>
            sym_grad_u_q = u_eval.get_symmetric_gradient(q);
            const VectorizedArray<number> div_u_q = trace(sym_grad_u_q);
            const VectorizedArray<number> val_p_q = p_eval.get_value(q);

            // Terms to be tested by phi_p:
            const VectorizedArray<number> pressure_terms =
              -cell_data->pressure_scaling * div_u_q;

            // Terms to be tested by the symmetric gradients of phi_u:
            SymmetricTensor<2,dim,VectorizedArray<number>>
            velocity_terms = viscosity_x_2 * sym_grad_u_q;

            for (unsigned int d=0; d<dim; ++d)
              velocity_terms[d][d] -= cell_data->pressure_scaling * val_p_q;

            if (cell_data->is_compressible)
              for (unsigned int d=0; d<dim; ++d)
                velocity_terms[d][d] -= viscosity_x_2 / 3. * div_u_q;

            // Add the Newton derivatives if required.
            if (cell_data->enable_newton_derivatives)
              {
                VectorizedArray<number> deta_deps_times_sym_grad_u(0.);
                VectorizedArray<number> eps_times_sym_grad_u(0.);
                VectorizedArray<number> deta_dp_times_p(0.);
                for (const unsigned int r : u_eval.quadrature_point_indices())
                  {
                    deta_deps_times_sym_grad_u += cell_data->newton_factor_wrt_strain_rate_table(cell,r)
                                                  * sym_grad_u[r];
                    deta_dp_times_p += cell_data->newton_factor_wrt_pressure_table(cell,r) * val_p[r];
                    if (cell_data->symmetrize_newton_system)
                      eps_times_sym_grad_u += cell_data->strain_rate_table(cell,r) * sym_grad_u[r];
                  }

                velocity_terms +=
                  ( cell_data->symmetrize_newton_system ?
                    ( cell_data->strain_rate_table(cell,q) * deta_deps_times_sym_grad_u +
                      cell_data->newton_factor_wrt_strain_rate_table(cell,q) * eps_times_sym_grad_u ) :
                    2. * cell_data->strain_rate_table(cell,q) * deta_deps_times_sym_grad_u )
                  +
                  2. * cell_data->strain_rate_table(cell,q) * deta_dp_times_p;
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

        if (cell_data->is_compressible)
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



  template <int dim>
  void StokesMatrixFreeHandler<dim>::declare_parameters(ParameterHandler &prm)
  {
    StokesMatrixFreeHandlerImplementation<dim,2>::declare_parameters(prm);
  }



  template <int dim, int velocity_degree>
  void
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::declare_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Solver parameters");
    prm.enter_subsection ("Matrix Free");
    {
      prm.declare_entry ("Output details", "false",
                         Patterns::Bool(),
                         "Turns on extra information for the matrix free GMG solver to be printed.");
      prm.declare_entry ("Execute solver timings", "false",
                         Patterns::Bool(),
                         "Executes different parts of the Stokes solver repeatedly and print timing information. "
                         "This is for internal benchmarking purposes: It is useful if you want to see how the solver "
                         "performs. Otherwise, you don't want to enable this, since it adds additional computational cost "
                         "to get the timing information.");
    }
    prm.leave_subsection ();
    prm.leave_subsection ();

  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim,velocity_degree>::parse_parameters(ParameterHandler &prm)
  {
    prm.enter_subsection ("Solver parameters");
    prm.enter_subsection ("Matrix Free");
    {
      print_details = prm.get_bool ("Output details");
      do_timings = prm.get_bool ("Execute solver timings");
    }
    prm.leave_subsection ();
    prm.leave_subsection ();
  }



  template <int dim, int velocity_degree>
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::StokesMatrixFreeHandlerImplementation (Simulator<dim> &simulator,
      ParameterHandler &prm)
    : sim(simulator),

      dof_handler_v(simulator.triangulation),
      dof_handler_p(simulator.triangulation),
      dof_handler_projection(simulator.triangulation),

      fe_v (FE_Q<dim>(sim.parameters.stokes_velocity_degree), dim),
      fe_p (FE_Q<dim>(sim.parameters.stokes_velocity_degree-1),1),

      // The finite element used to describe the viscosity on the active level
      // and to project the viscosity to GMG levels needs to be DGQ1 if we are
      // using a degree 1 representation of viscosity, and DGQ0 if we are using
      // a cellwise constant average.
      fe_projection(FE_DGQ<dim>(sim.parameters.material_averaging
                                ==
                                MaterialModel::MaterialAveraging::AveragingOperation::project_to_Q1
                                ||
                                sim.parameters.material_averaging
                                ==
                                MaterialModel::MaterialAveraging::AveragingOperation::project_to_Q1_only_viscosity
                                ? 1 : 0), 1)
  {
    parse_parameters(prm);
    CitationInfo::add("mf");

    // Sorry, not any time soon:
    AssertThrow(!sim.parameters.include_melt_transport, ExcNotImplemented());
    // Not very difficult to do, but will require a different mass matrix
    // operator:
    AssertThrow(!sim.parameters.use_locally_conservative_discretization, ExcNotImplemented());

    // sanity check:
    Assert(sim.introspection.variable("velocity").block_index==0, ExcNotImplemented());
    Assert(sim.introspection.variable("pressure").block_index==1, ExcNotImplemented());

    // We currently only support averaging of the viscosity to a constant or Q1:
    using avg = MaterialModel::MaterialAveraging::AveragingOperation;
    AssertThrow((sim.parameters.material_averaging &
                 (avg::arithmetic_average | avg::harmonic_average | avg::geometric_average
                  | avg::pick_largest | avg::project_to_Q1 | avg::log_average
                  | avg::harmonic_average_only_viscosity | avg::geometric_average_only_viscosity
                  | avg::project_to_Q1_only_viscosity)) != 0,
                ExcMessage("The matrix-free Stokes solver currently only works if material model averaging "
                           "is enabled. If no averaging is desired, consider using ``project to Q1 only "
                           "viscosity''."));

    // Currently cannot solve compressible flow with implicit reference density
    if (sim.material_model->is_compressible() == true)
      AssertThrow(sim.parameters.formulation_mass_conservation !=
                  Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile,
                  ExcNotImplemented());
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::assemble ()
  {
    if (sim.mesh_deformation)
      {
        // Update the geometry information stored in the MatrixFree
        // objects from the mapping.  Grab the mapping stored in the
        // object and do not replace with sim.mapping as we have
        // different mappings per level.
        for (auto &obj : matrix_free_objects)
          obj->update_mapping(*obj->get_mapping_info().mapping);
      }

    evaluate_material_model();

    correct_stokes_rhs();
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::evaluate_material_model ()
  {
    dealii::LinearAlgebra::distributed::Vector<double> active_viscosity_vector(dof_handler_projection.locally_owned_dofs(),
                                                                               sim.triangulation.get_communicator());

    const Quadrature<dim> &quadrature_formula = sim.introspection.quadratures.velocities;

    double minimum_viscosity_local = std::numeric_limits<double>::max();
    double maximum_viscosity_local = std::numeric_limits<double>::lowest();

    // Fill the DGQ0 or DGQ1 vector of viscosity values on the active mesh
    {
      FEValues<dim> fe_values (*sim.mapping,
                               sim.finite_element,
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, sim.introspection.n_compositional_fields);
      in.requested_properties = MaterialModel::MaterialProperties::viscosity;
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, sim.introspection.n_compositional_fields);

      // This function call computes a cellwise projection of data defined at quadrature points to
      // a vector defined by the projection DoFHandler. As an input, we must define a lambda which returns
      // a viscosity value for each quadrature point of the given cell. The projection is then stored in
      // the active level viscosity vector provided.
      Utilities::project_cellwise<dim, dealii::LinearAlgebra::distributed::Vector<double>>(*(sim.mapping),
          dof_handler_projection,
          0,
          quadrature_formula,
          [&](const typename DoFHandler<dim>::active_cell_iterator & cell,
              const std::vector<Point<dim>> & /*q_points*/,
              std::vector<double> &values) -> void
      {
        typename DoFHandler<dim>::active_cell_iterator FEQ_cell(&sim.triangulation,
        cell->level(),
        cell->index(),
        &(sim.dof_handler));

        fe_values.reinit (FEQ_cell);
        in.reinit(fe_values, FEQ_cell, sim.introspection, sim.current_linearization_point);

        // Query the material model for the active level viscosities
        sim.material_model->fill_additional_material_model_inputs(in, sim.current_linearization_point, fe_values, sim.introspection);
        sim.material_model->evaluate(in, out);

        // If using a cellwise average for viscosity, average the values here.
        // When the projection is computed, this will set the viscosity exactly
        // to this averaged value.
        if (dof_handler_projection.get_fe().degree == 0)
          MaterialModel::MaterialAveraging::average (sim.parameters.material_averaging,
          FEQ_cell,
          quadrature_formula,
          *sim.mapping,
          in.requested_properties,
          out);

        for (unsigned int i=0; i<values.size(); ++i)
          {
            // Find the local max/min of the evaluated viscosities.
            minimum_viscosity_local = std::min(minimum_viscosity_local, out.viscosities[i]);
            maximum_viscosity_local = std::max(maximum_viscosity_local, out.viscosities[i]);

            values[i] = out.viscosities[i];
          }
        return;
      },
      active_viscosity_vector);

      active_viscosity_vector.compress(VectorOperation::insert);
    }

    minimum_viscosity = dealii::Utilities::MPI::min(minimum_viscosity_local, sim.triangulation.get_communicator());
    maximum_viscosity = dealii::Utilities::MPI::max(maximum_viscosity_local, sim.triangulation.get_communicator());

    FEValues<dim> fe_values_projection (*(sim.mapping),
                                        fe_projection,
                                        quadrature_formula,
                                        update_values);

    // Create active mesh viscosity table.
    {

      const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_cell_batches();

      const unsigned int n_q_points = quadrature_formula.size();

      std::vector<double> values_on_quad;

      // One value per cell is required for DGQ0 projection and n_q_points
      // values per cell for DGQ1.
      if (dof_handler_projection.get_fe().degree == 0)
        active_cell_data.viscosity.reinit(TableIndices<2>(n_cells, 1));
      else if (dof_handler_projection.get_fe().degree == 1)
        {
          values_on_quad.resize(n_q_points);
          active_cell_data.viscosity.reinit(TableIndices<2>(n_cells, n_q_points));
        }
      else
        Assert(false, ExcInternalError());

      std::vector<types::global_dof_index> local_dof_indices(fe_projection.dofs_per_cell);
      for (unsigned int cell=0; cell<n_cells; ++cell)
        {
          const unsigned int n_components_filled = stokes_matrix.get_matrix_free()->n_active_entries_per_cell_batch(cell);

          for (unsigned int i=0; i<n_components_filled; ++i)
            {
              typename DoFHandler<dim>::active_cell_iterator FEQ_cell =
                stokes_matrix.get_matrix_free()->get_cell_iterator(cell,i);
              typename DoFHandler<dim>::active_cell_iterator DG_cell(&(sim.triangulation),
                                                                     FEQ_cell->level(),
                                                                     FEQ_cell->index(),
                                                                     &dof_handler_projection);
              DG_cell->get_active_or_mg_dof_indices(local_dof_indices);

#ifdef DEBUG
              {
                // Verify that all MatrixFree objects iterate over cells in the same way:
                typename DoFHandler<dim>::active_cell_iterator s_cell =
                  Schur_complement_block_matrix.get_matrix_free()->get_cell_iterator(cell,i,1);
                double distance_s = s_cell->center().distance(FEQ_cell->center());
                Assert(distance_s < 1e-10, ExcInternalError());

                typename DoFHandler<dim>::active_cell_iterator A_cell =
                  A_block_matrix.get_matrix_free()->get_cell_iterator(cell,i);
                double distance_A = A_cell->center().distance(FEQ_cell->center());
                Assert(distance_A < 1e-10, ExcInternalError());
              }
#endif

              // For DGQ0, we simply use the viscosity at the single
              // support point of the element. For DGQ1, we must project
              // back to quadrature point values.
              if (dof_handler_projection.get_fe().degree == 0)
                active_cell_data.viscosity(cell, 0)[i] = active_viscosity_vector(local_dof_indices[0]);
              else
                {
                  fe_values_projection.reinit(DG_cell);
                  fe_values_projection.get_function_values(active_viscosity_vector,
                                                           local_dof_indices,
                                                           values_on_quad);

                  // Do not allow viscosity to be greater than or less than the limits
                  // of the evaluated viscosity on the active level.
                  for (unsigned int q=0; q<n_q_points; ++q)
                    active_cell_data.viscosity(cell, q)[i]
                      = std::min(std::max(values_on_quad[q], minimum_viscosity), maximum_viscosity);
                }
            }
        }
    }

    active_cell_data.is_compressible = sim.material_model->is_compressible();
    active_cell_data.pressure_scaling = sim.pressure_scaling;

    // Store viscosity tables and other data into the active level matrix-free objects.
    stokes_matrix.set_cell_data(active_cell_data);

    if (sim.parameters.n_expensive_stokes_solver_steps > 0)
      {
        A_block_matrix.set_cell_data(active_cell_data);
        Schur_complement_block_matrix.set_cell_data(active_cell_data);
      }

    const unsigned int n_levels = sim.triangulation.n_global_levels();
    level_cell_data.resize(0,n_levels-1);

    MGLevelObject<dealii::LinearAlgebra::distributed::Vector<GMGNumberType>> level_viscosity_vector;
    level_viscosity_vector.resize(0,n_levels-1);

    // Project the active level viscosity vector to multilevel vector representations
    // using MG transfer objects. This transfer is based on the same linear operator used to
    // transfer data inside a v-cycle.
    MGTransferMF<dim,GMGNumberType> transfer;

    transfer.build(dof_handler_projection);

    transfer.interpolate_to_mg(dof_handler_projection,
                               level_viscosity_vector,
                               active_viscosity_vector);

    for (unsigned int level=0; level<n_levels; ++level)
      {
        level_cell_data[level].is_compressible = sim.material_model->is_compressible();
        level_cell_data[level].pressure_scaling = sim.pressure_scaling;

        // Create viscosity tables on each level.
        const unsigned int n_cells = mg_matrices_A_block[level].get_matrix_free()->n_cell_batches();

        const unsigned int n_q_points = quadrature_formula.size();

        std::vector<GMGNumberType> values_on_quad;

        // One value per cell is required for DGQ0 projection and n_q_points
        // values per cell for DGQ1.
        if (dof_handler_projection.get_fe().degree == 0)
          level_cell_data[level].viscosity.reinit(TableIndices<2>(n_cells, 1));
        else
          {
            values_on_quad.resize(n_q_points);
            level_cell_data[level].viscosity.reinit(TableIndices<2>(n_cells, n_q_points));
          }

        std::vector<types::global_dof_index> local_dof_indices(fe_projection.dofs_per_cell);
        for (unsigned int cell=0; cell<n_cells; ++cell)
          {
            const unsigned int n_components_filled = mg_matrices_A_block[level].get_matrix_free()->n_active_entries_per_cell_batch(cell);

            for (unsigned int i=0; i<n_components_filled; ++i)
              {
                typename DoFHandler<dim>::level_cell_iterator FEQ_cell =
                  mg_matrices_A_block[level].get_matrix_free()->get_cell_iterator(cell,i);
                typename DoFHandler<dim>::level_cell_iterator DG_cell(&(sim.triangulation),
                                                                      FEQ_cell->level(),
                                                                      FEQ_cell->index(),
                                                                      &dof_handler_projection);
                DG_cell->get_active_or_mg_dof_indices(local_dof_indices);

                // For DGQ0, we simply use the viscosity at the single
                // support point of the element. For DGQ1, we must project
                // back to quadrature point values.
                if (dof_handler_projection.get_fe().degree == 0)
                  level_cell_data[level].viscosity(cell, 0)[i] = level_viscosity_vector[level](local_dof_indices[0]);
                else
                  {
                    fe_values_projection.reinit(DG_cell);
                    fe_values_projection.get_function_values(level_viscosity_vector[level],
                                                             local_dof_indices,
                                                             values_on_quad);

                    // Do not allow viscosity to be greater than or less than the limits
                    // of the evaluated viscosity on the active level.
                    for (unsigned int q=0; q<n_q_points; ++q)
                      level_cell_data[level].viscosity(cell,q)[i]
                        = std::min(std::max(values_on_quad[q], static_cast<GMGNumberType>(minimum_viscosity)),
                                   static_cast<GMGNumberType>(maximum_viscosity));
                  }
              }
          }

        // Store viscosity tables and other data into the multigrid level matrix-free objects.
        mg_matrices_A_block[level].set_cell_data (level_cell_data[level]);
        mg_matrices_Schur_complement[level].set_cell_data (level_cell_data[level]);
      }

    {
      // create active mesh tables for derivatives needed in Newton method
      // and the strain rate.
      if (sim.newton_handler != nullptr
          && sim.newton_handler->parameters.newton_derivative_scaling_factor != 0)
        {
          const double newton_derivative_scaling_factor =
            sim.newton_handler->parameters.newton_derivative_scaling_factor;

          active_cell_data.enable_newton_derivatives = true;

          // TODO: these are not implemented yet
          for (unsigned int level=0; level<n_levels; ++level)
            level_cell_data[level].enable_newton_derivatives = false;


          FEValues<dim> fe_values (*sim.mapping,
                                   sim.finite_element,
                                   quadrature_formula,
                                   update_values   |
                                   update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

          MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, sim.introspection.n_compositional_fields);
          MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, sim.introspection.n_compositional_fields);
          sim.newton_handler->create_material_model_outputs(out);
          if (sim.parameters.enable_elasticity &&
              out.template get_additional_output<MaterialModel::ElasticOutputs<dim>>() == nullptr)
            out.additional_outputs.push_back(std::make_unique<MaterialModel::ElasticOutputs<dim>>(out.n_evaluation_points()));

          const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_cell_batches();
          const unsigned int n_q_points = quadrature_formula.size();

          active_cell_data.strain_rate_table.reinit(TableIndices<2>(n_cells, n_q_points));
          active_cell_data.newton_factor_wrt_pressure_table.reinit(TableIndices<2>(n_cells, n_q_points));
          active_cell_data.newton_factor_wrt_strain_rate_table.reinit(TableIndices<2>(n_cells, n_q_points));

          for (unsigned int cell=0; cell<n_cells; ++cell)
            {
              const unsigned int n_components_filled = stokes_matrix.get_matrix_free()->n_active_entries_per_cell_batch(cell);

              for (unsigned int i=0; i<n_components_filled; ++i)
                {
                  typename DoFHandler<dim>::active_cell_iterator matrix_free_cell =
                    stokes_matrix.get_matrix_free()->get_cell_iterator(cell,i);
                  typename DoFHandler<dim>::active_cell_iterator simulator_cell(&(sim.triangulation),
                                                                                matrix_free_cell->level(),
                                                                                matrix_free_cell->index(),
                                                                                &(sim.dof_handler));

                  fe_values.reinit(simulator_cell);
                  in.reinit(fe_values, simulator_cell, sim.introspection, sim.current_linearization_point);

                  sim.material_model->fill_additional_material_model_inputs(in, sim.current_linearization_point, fe_values, sim.introspection);
                  sim.material_model->evaluate(in, out);

                  MaterialModel::MaterialAveraging::average(sim.parameters.material_averaging,
                                                            in.current_cell,
                                                            fe_values.get_quadrature(),
                                                            *sim.mapping,
                                                            in.requested_properties,
                                                            out);

                  Assert(std::isfinite(in.strain_rate[0].norm()),
                         ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                    "not filled by the caller."));

                  const MaterialModel::MaterialModelDerivatives<dim> *derivatives
                    = out.template get_additional_output<MaterialModel::MaterialModelDerivatives<dim>>();

                  Assert(derivatives != nullptr,
                         ExcMessage ("Error: The Newton method requires the material to "
                                     "compute derivatives."));

                  const MaterialModel::ElasticOutputs<dim> *elastic_out
                    = out.template get_additional_output<MaterialModel::ElasticOutputs<dim>>();

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      // use the correct strain rate for the Jacobian
                      // when elasticity is enabled use viscoelastic strain rate
                      // when stabilization is enabled, use the deviatoric strain rate because the SPD factor
                      // that is computed is only safe for the deviatoric strain rate (see PR #5580 and issue #5555)
                      SymmetricTensor<2,dim> effective_strain_rate = in.strain_rate[q];
                      if (elastic_out != nullptr)
                        effective_strain_rate = elastic_out->viscoelastic_strain_rate[q];
                      else if ((sim.newton_handler->parameters.velocity_block_stabilization & Newton::Parameters::Stabilization::PD) != Newton::Parameters::Stabilization::none)
                        effective_strain_rate = deviator(effective_strain_rate);

                      // use the spd factor when the stabilization is PD or SPD.
                      const double alpha =  (sim.newton_handler->parameters.velocity_block_stabilization
                                             & Newton::Parameters::Stabilization::PD)
                                            != Newton::Parameters::Stabilization::none
                                            ?
                                            Utilities::compute_spd_factor<dim>(out.viscosities[q],
                                                                               effective_strain_rate,
                                                                               derivatives->viscosity_derivative_wrt_strain_rate[q],
                                                                               sim.newton_handler->parameters.SPD_safety_factor)
                                            :
                                            1.0;

                      active_cell_data.newton_factor_wrt_pressure_table(cell,q)[i]
                        = derivatives->viscosity_derivative_wrt_pressure[q] *
                          derivatives->viscosity_derivative_averaging_weights[q] *
                          newton_derivative_scaling_factor;
                      Assert(std::isfinite(active_cell_data.newton_factor_wrt_pressure_table(cell,q)[i]),
                             ExcMessage("active_cell_data.newton_factor_wrt_pressure_table is not finite: " + std::to_string(active_cell_data.newton_factor_wrt_pressure_table(cell,q)[i]) +
                                        ". Relevant variables are derivatives->viscosity_derivative_wrt_pressure[q] = " + std::to_string(derivatives->viscosity_derivative_wrt_pressure[q]) +
                                        ", derivatives->viscosity_derivative_averaging_weights[q] = " + std::to_string(derivatives->viscosity_derivative_averaging_weights[q]) +
                                        ", and newton_derivative_scaling_factor = " + std::to_string(newton_derivative_scaling_factor)));

                      for (unsigned int m=0; m<dim; ++m)
                        for (unsigned int n=0; n<dim; ++n)
                          {
                            active_cell_data.strain_rate_table(cell, q)[m][n][i]
                              = effective_strain_rate[m][n];

                            active_cell_data.newton_factor_wrt_strain_rate_table(cell, q)[m][n][i]
                              = derivatives->viscosity_derivative_wrt_strain_rate[q][m][n] *
                                derivatives->viscosity_derivative_averaging_weights[q] *
                                newton_derivative_scaling_factor * alpha;

                            Assert(std::isfinite(active_cell_data.strain_rate_table(cell, q)[m][n][i]),
                                   ExcMessage("active_cell_data.strain_rate_table has an element which is not finite: " + std::to_string(active_cell_data.strain_rate_table(cell, q)[m][n][i])));
                            Assert(std::isfinite(active_cell_data.newton_factor_wrt_strain_rate_table(cell, q)[m][n][i]),
                                   ExcMessage("active_cell_data.newton_factor_wrt_strain_rate_table has an element which is not finite: " + std::to_string(active_cell_data.newton_factor_wrt_strain_rate_table(cell, q)[m][n][i])));
                          }
                    }
                }
            }

          // symmetrize the Newton_system when the stabilization is symmetric or SPD
          const bool symmetrize_newton_system =
            (sim.newton_handler->parameters.velocity_block_stabilization & Newton::Parameters::Stabilization::symmetric)
            != Newton::Parameters::Stabilization::none;
          active_cell_data.symmetrize_newton_system = symmetrize_newton_system;
        }
      else
        {
          // delete data used for Newton derivatives if necessary
          // TODO: use Table::clear() once implemented in 10.0.pre
          active_cell_data.enable_newton_derivatives = false;
          active_cell_data.newton_factor_wrt_pressure_table.reinit(TableIndices<2>(0,0));
          active_cell_data.strain_rate_table.reinit(TableIndices<2>(0,0));
          active_cell_data.newton_factor_wrt_strain_rate_table.reinit(TableIndices<2>(0,0));

          for (unsigned int level=0; level<n_levels; ++level)
            level_cell_data[level].enable_newton_derivatives = false;
        }
    }

    {
      // Create active mesh tables to store the product of the pressure perturbation and
      // the normalized gravity used in the free surface stabilization.
      // Currently, mutilevel is not implemented yet, it may slow down the convergence.

      // TODO: implement multilevel surface terms for the free surface stabilization.

      active_cell_data.apply_stabilization_free_surface_faces = sim.mesh_deformation
                                                                && !sim.mesh_deformation->get_free_surface_boundary_indicators().empty();
      if (active_cell_data.apply_stabilization_free_surface_faces == true)
        {
          const double free_surface_theta = sim.mesh_deformation->get_free_surface_theta();

          const Quadrature<dim-1> &face_quadrature_formula = sim.introspection.face_quadratures.velocities;

          const unsigned int n_face_q_points = face_quadrature_formula.size();

          // We need the gradients for the material model inputs.
          FEFaceValues<dim> fe_face_values (*sim.mapping,
                                            sim.finite_element,
                                            face_quadrature_formula,
                                            update_values   |
                                            update_gradients |
                                            update_quadrature_points |
                                            update_JxW_values);

          const unsigned int n_faces_boundary = stokes_matrix.get_matrix_free()->n_boundary_face_batches();
          const unsigned int n_faces_interior = stokes_matrix.get_matrix_free()->n_inner_face_batches();

          active_cell_data.free_surface_boundary_indicators =
            sim.mesh_deformation->get_free_surface_boundary_indicators();

          MaterialModel::MaterialModelInputs<dim> face_material_inputs(n_face_q_points, sim.introspection.n_compositional_fields);
          face_material_inputs.requested_properties = MaterialModel::MaterialProperties::density;
          MaterialModel::MaterialModelOutputs<dim> face_material_outputs(n_face_q_points, sim.introspection.n_compositional_fields);

          active_cell_data.free_surface_stabilization_term_table.reinit(n_faces_boundary, n_face_q_points);

          for (unsigned int face=n_faces_interior; face<n_faces_boundary + n_faces_interior; ++face)
            {
              const unsigned int n_components_filled = stokes_matrix.get_matrix_free()->n_active_entries_per_face_batch(face);

              for (unsigned int i=0; i<n_components_filled; ++i)
                {
                  // The first element of the pair is the active cell iterator
                  // the second element of the pair is the face number
                  const auto cell_face_pair = stokes_matrix.get_matrix_free()->get_face_iterator(face, i, true);

                  typename DoFHandler<dim>::active_cell_iterator matrix_free_cell =
                    cell_face_pair.first;
                  typename DoFHandler<dim>::active_cell_iterator simulator_cell(&(sim.triangulation),
                                                                                matrix_free_cell->level(),
                                                                                matrix_free_cell->index(),
                                                                                &(sim.dof_handler));

                  const types::boundary_id boundary_indicator = stokes_matrix.get_matrix_free()->get_boundary_id(face);
                  Assert(boundary_indicator == simulator_cell->face(cell_face_pair.second)->boundary_id(), ExcInternalError());

                  // only apply on free surface faces
                  if (active_cell_data.free_surface_boundary_indicators.find(boundary_indicator)
                      == active_cell_data.free_surface_boundary_indicators.end())
                    continue;

                  fe_face_values.reinit(simulator_cell, cell_face_pair.second);

                  face_material_inputs.reinit  (fe_face_values,
                                                simulator_cell,
                                                sim.introspection,
                                                sim.solution);
                  face_material_inputs.requested_properties = MaterialModel::MaterialProperties::density;
                  sim.material_model->evaluate(face_material_inputs, face_material_outputs);

                  for (unsigned int q = 0; q < n_face_q_points; ++q)
                    {
                      const Tensor<1,dim>
                      gravity = sim.gravity_model->gravity_vector(fe_face_values.quadrature_point(q));
                      const double g_norm = gravity.norm();

                      const Tensor<1,dim> g_hat = (g_norm == 0.0 ? Tensor<1,dim>() : gravity/g_norm);

                      const double pressure_perturbation = face_material_outputs.densities[q] *
                                                           sim.time_step *
                                                           free_surface_theta *
                                                           g_norm;
                      for (unsigned int d = 0; d < dim; ++d)
                        active_cell_data.free_surface_stabilization_term_table(face - n_faces_interior, q)[d][i]
                          = pressure_perturbation * g_hat[d];
                    }
                }
            }
        }
    }
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::correct_stokes_rhs()
  {
    // We never include Newton terms in step 0 and after that we solve with zero boundary conditions.
    // Therefore, we don't need to include Newton terms here.

    const bool is_compressible = sim.material_model->is_compressible();

    dealii::LinearAlgebra::distributed::BlockVector<double> rhs_correction(2);
    dealii::LinearAlgebra::distributed::BlockVector<double> u0(2);

    stokes_matrix.initialize_dof_vector(rhs_correction);
    stokes_matrix.initialize_dof_vector(u0);

    // The vector u0 is a zero vector, but we need to ensure that it
    // has the correct boundary values:
    u0 = 0;

#if DEAL_II_VERSION_GTE(9,6,0)
    IndexSet stokes_dofs (sim.dof_handler.n_dofs());
    stokes_dofs.add_range (0, u0.size());
    const AffineConstraints<double> current_stokes_constraints
      = sim.current_constraints.get_view (stokes_dofs);
    current_stokes_constraints.distribute(u0);
#else
    sim.current_constraints.distribute(u0);
#endif

    u0.update_ghost_values();

    rhs_correction = 0;

    FEEvaluation<dim,velocity_degree,velocity_degree+1,dim,double>
    velocity (*stokes_matrix.get_matrix_free(), 0);
    FEEvaluation<dim,velocity_degree-1,velocity_degree+1,1,double>
    pressure (*stokes_matrix.get_matrix_free(), 1);

    const bool use_viscosity_at_quadrature_points
      = (active_cell_data.viscosity.size(1) == velocity.n_q_points);

    const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_cell_batches();

    // Much like the matrix-free apply_add() functions compute a matrix-vector
    // product by looping over cells and applying local matrix operations,
    // here we apply the negative of the stokes_matrix operator to u0.
    for (unsigned int cell=0; cell<n_cells; ++cell)
      {
        VectorizedArray<double> viscosity_x_2 = 2.0*active_cell_data.viscosity(cell, 0);

        // We must use read_dof_values_plain() as to not overwrite boundary information
        // with the zero boundary used by the stokes_matrix operator.
        velocity.reinit (cell);
        velocity.read_dof_values_plain (u0.block(0));
        velocity.evaluate (EvaluationFlags::gradients);

        pressure.reinit (cell);
        pressure.read_dof_values_plain (u0.block(1));
        pressure.evaluate (EvaluationFlags::values);

        for (const unsigned int q : velocity.quadrature_point_indices())
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              viscosity_x_2 = 2.0*active_cell_data.viscosity(cell, q);

            SymmetricTensor<2,dim,VectorizedArray<double>> sym_grad_u =
              velocity.get_symmetric_gradient (q);
            const VectorizedArray<double> pres = pressure.get_value(q);
            const VectorizedArray<double> div = trace(sym_grad_u);
            pressure.submit_value(sim.pressure_scaling*div, q);

            sym_grad_u *= viscosity_x_2;

            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= sim.pressure_scaling*pres;

            if (is_compressible)
              for (unsigned int d=0; d<dim; ++d)
                sym_grad_u[d][d] -= viscosity_x_2/3.0*div;

            velocity.submit_symmetric_gradient(-1.0*sym_grad_u, q);
          }


        velocity.integrate_scatter (EvaluationFlags::gradients,
                                    rhs_correction.block(0));

        pressure.integrate_scatter (EvaluationFlags::values,
                                    rhs_correction.block(1));
      }

    if (active_cell_data.apply_stabilization_free_surface_faces)
      {
        const unsigned int n_faces_boundary = stokes_matrix.get_matrix_free()->n_boundary_face_batches();
        const unsigned int n_faces_interior = stokes_matrix.get_matrix_free()->n_inner_face_batches();

        FEFaceEvaluation<dim,velocity_degree,velocity_degree+1,dim,double>
        velocity_boundary(*stokes_matrix.get_matrix_free());

        for (unsigned int face=n_faces_interior; face<n_faces_boundary + n_faces_interior; ++face)
          {
            const auto boundary_id = stokes_matrix.get_matrix_free()->get_boundary_id(face);
            if (active_cell_data.free_surface_boundary_indicators.find(boundary_id)
                == active_cell_data.free_surface_boundary_indicators.end())
              continue;

            velocity_boundary.reinit(face);
            velocity_boundary.read_dof_values_plain (u0.block(0));
            velocity_boundary.evaluate (EvaluationFlags::values);

            for (const unsigned int q : velocity_boundary.quadrature_point_indices())
              {
                const Tensor<1, dim, VectorizedArray<double>> phi_u_i = velocity_boundary.get_value(q);
#if DEAL_II_VERSION_GTE(9, 7, 0)
                const auto &normal_vector = velocity_boundary.normal_vector(q);
#else
                const auto &normal_vector = velocity_boundary.get_normal_vector(q);
#endif
                const auto stabilization_tensor = active_cell_data.free_surface_stabilization_term_table(face - n_faces_interior, q);
                const auto value_submit = (stabilization_tensor * phi_u_i) * normal_vector;
                velocity_boundary.submit_value(value_submit, q);

              }
            velocity_boundary.integrate_scatter(EvaluationFlags::values,
                                                rhs_correction.block(0));
          }
      }

    rhs_correction.compress(VectorOperation::add);

    // Copy to the correct vector type and add the correction to the system rhs.
    LinearAlgebra::BlockVector stokes_rhs_correction (sim.introspection.index_sets.stokes_partitioning, sim.mpi_communicator);
    internal::ChangeVectorTypes::copy(stokes_rhs_correction,rhs_correction);

    sim.system_rhs.block(0) += stokes_rhs_correction.block(0);
    sim.system_rhs.block(1) += stokes_rhs_correction.block(1);
  }



  template <int dim, int velocity_degree>
  std::pair<double,double> StokesMatrixFreeHandlerImplementation<dim,velocity_degree>::solve(LinearAlgebra::BlockVector &solution_vector)
  {
    double initial_nonlinear_residual = numbers::signaling_nan<double>();
    double final_linear_residual      = numbers::signaling_nan<double>();

    // Below we define all the objects needed to build the GMG preconditioner:
    using VectorType = dealii::LinearAlgebra::distributed::Vector<GMGNumberType>;

    // ABlock GMG Smoother: Chebyshev, degree 4. Parameter values were chosen
    // by trial and error. We use a more powerful version of the smoother on the
    // coarsest level than on the other levels.
    using ASmootherType = PreconditionChebyshev<GMGABlockMatrixType,VectorType>;
    mg::SmootherRelaxation<ASmootherType, VectorType>
    mg_smoother_A;
    {
      MGLevelObject<typename ASmootherType::AdditionalData> smoother_data_A;
      smoother_data_A.resize(0, sim.triangulation.n_global_levels()-1);
      for (unsigned int level = 0; level<sim.triangulation.n_global_levels(); ++level)
        {
          if (level > 0)
            {
              smoother_data_A[level].smoothing_range = 15.;
              smoother_data_A[level].degree = 4;
              smoother_data_A[level].eig_cg_n_iterations = 10;
            }
          else
            {
              smoother_data_A[0].smoothing_range = 1e-3;
              smoother_data_A[0].degree = 8;
              smoother_data_A[0].eig_cg_n_iterations = 100;
            }
          smoother_data_A[level].preconditioner = mg_matrices_A_block[level].get_matrix_diagonal_inverse();
        }
      mg_smoother_A.initialize(mg_matrices_A_block, smoother_data_A);
    }

    // Schur complement matrix GMG Smoother: Chebyshev, degree 4. Parameter values
    // were chosen by trial and error. We use a more powerful version of the smoother
    // on the coarsest level than on the other levels.
    using MSmootherType = PreconditionChebyshev<GMGSchurComplementMatrixType,VectorType>;
    mg::SmootherRelaxation<MSmootherType, VectorType>
    mg_smoother_Schur(4);
    {
      MGLevelObject<typename MSmootherType::AdditionalData> smoother_data_Schur;
      smoother_data_Schur.resize(0, sim.triangulation.n_global_levels()-1);
      for (unsigned int level = 0; level<sim.triangulation.n_global_levels(); ++level)
        {
          if (level > 0)
            {
              smoother_data_Schur[level].smoothing_range = 15.;
              smoother_data_Schur[level].degree = 4;
              smoother_data_Schur[level].eig_cg_n_iterations = 10;
            }
          else
            {
              smoother_data_Schur[0].smoothing_range = 1e-3;
              smoother_data_Schur[0].degree = 8;
              smoother_data_Schur[0].eig_cg_n_iterations = 100;
            }
          smoother_data_Schur[level].preconditioner = mg_matrices_Schur_complement[level].get_matrix_diagonal_inverse();
        }
      mg_smoother_Schur.initialize(mg_matrices_Schur_complement, smoother_data_Schur);
    }

    // Estimate the eigenvalues for the Chebyshev smoothers.

    types::global_dof_index coarse_A_size = numbers::invalid_dof_index, coarse_S_size = numbers::invalid_dof_index;

    //TODO: The setup for the smoother (as well as the entire GMG setup) should
    //       be moved to an assembly timing block instead of the Stokes solve
    //       timing block (as is currently the case).
    for (unsigned int level = 0; level<sim.triangulation.n_global_levels(); ++level)
      {
        VectorType temp_velocity;
        VectorType temp_pressure;
        mg_matrices_A_block[level].initialize_dof_vector(temp_velocity);
        mg_matrices_Schur_complement[level].initialize_dof_vector(temp_pressure);

        mg_smoother_A[level].estimate_eigenvalues(temp_velocity);
        mg_smoother_Schur[level].estimate_eigenvalues(temp_pressure);

        if (level==0)
          {
            coarse_A_size = temp_velocity.size();
            coarse_S_size = temp_pressure.size();
          }
      }


    // Coarse Solver is just an application of the Chebyshev smoother setup
    // in such a way to be a solver
    //ABlock GMG
    MGCoarseGridApplySmoother<VectorType> mg_coarse_A;
    mg_coarse_A.initialize(mg_smoother_A);

    //Schur complement matrix GMG
    MGCoarseGridApplySmoother<VectorType> mg_coarse_Schur;
    mg_coarse_Schur.initialize(mg_smoother_Schur);


    if (print_details)
      {
        sim.pcout << std::endl
                  << "    GMG coarse size A: " << coarse_A_size << ", coarse size S: " << coarse_S_size << std::endl
                  << "    GMG n_levels: " << sim.triangulation.n_global_levels() << std::endl
                  << "    Viscosity range: " << minimum_viscosity << " - " << maximum_viscosity << std::endl;

        const double imbalance = MGTools::workload_imbalance(sim.triangulation);
        sim.pcout << "    GMG workload imbalance: " << imbalance << std::endl
                  << "    Stokes solver: " << std::flush;
      }

    // Interface matrices
    // Ablock GMG
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<GMGABlockMatrixType>> mg_interface_matrices_A;
    mg_interface_matrices_A.resize(0, sim.triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<sim.triangulation.n_global_levels(); ++level)
      mg_interface_matrices_A[level].initialize(mg_matrices_A_block[level]);
    mg::Matrix<VectorType> mg_interface_A(mg_interface_matrices_A);

    // Schur complement matrix GMG
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<GMGSchurComplementMatrixType>> mg_interface_matrices_Schur;
    mg_interface_matrices_Schur.resize(0, sim.triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<sim.triangulation.n_global_levels(); ++level)
      mg_interface_matrices_Schur[level].initialize(mg_matrices_Schur_complement[level]);
    mg::Matrix<VectorType> mg_interface_Schur(mg_interface_matrices_Schur);

    // MG Matrix
    mg::Matrix<VectorType> mg_matrix_A(mg_matrices_A_block);
    mg::Matrix<VectorType> mg_matrix_Schur(mg_matrices_Schur_complement);

    // MG object
    // ABlock GMG
    Multigrid<VectorType> mg_A(mg_matrix_A,
                               mg_coarse_A,
                               mg_transfer_A_block,
                               mg_smoother_A,
                               mg_smoother_A);
    mg_A.set_edge_matrices(mg_interface_A, mg_interface_A);

    // Schur complement matrix GMG
    Multigrid<VectorType> mg_Schur(mg_matrix_Schur,
                                   mg_coarse_Schur,
                                   mg_transfer_Schur_complement,
                                   mg_smoother_Schur,
                                   mg_smoother_Schur);
    mg_Schur.set_edge_matrices(mg_interface_Schur, mg_interface_Schur);

    // GMG Preconditioner for ABlock and Schur complement
    using GMGPreconditioner = PreconditionMG<dim, VectorType, MGTransferMF<dim,GMGNumberType>>;
    GMGPreconditioner prec_A(dof_handler_v, mg_A, mg_transfer_A_block);
    GMGPreconditioner prec_Schur(dof_handler_p, mg_Schur, mg_transfer_Schur_complement);


    // Many parts of the solver depend on the block layout (velocity = 0,
    // pressure = 1). For example the linearized_stokes_initial_guess vector or the StokesBlock matrix
    // wrapper. Let us make sure that this holds (and shorten their names):
    const unsigned int block_vel = sim.introspection.block_indices.velocities;
    const unsigned int block_p = (sim.parameters.include_melt_transport) ?
                                 sim.introspection.variable("fluid pressure").block_index
                                 : sim.introspection.block_indices.pressure;

    LinearAlgebra::BlockVector distributed_stokes_solution (sim.introspection.index_sets.stokes_partitioning,
                                                            sim.mpi_communicator);
    // extract Stokes parts of rhs vector
    LinearAlgebra::BlockVector distributed_stokes_rhs(sim.introspection.index_sets.stokes_partitioning,
                                                      sim.mpi_communicator);

    distributed_stokes_rhs.block(block_vel) = sim.system_rhs.block(block_vel);
    distributed_stokes_rhs.block(block_p) = sim.system_rhs.block(block_p);

    Assert(block_vel == 0, ExcNotImplemented());
    Assert(block_p == 1, ExcNotImplemented());
    Assert(!sim.parameters.include_melt_transport
           || sim.introspection.variable("compaction pressure").block_index == 1, ExcNotImplemented());

    // create a completely distributed vector that will be used for
    // the scaled and denormalized solution and later used as a
    // starting guess for the linear solver
    LinearAlgebra::BlockVector linearized_stokes_initial_guess (sim.introspection.index_sets.stokes_partitioning,
                                                                sim.mpi_communicator);

    // copy the velocity and pressure from current_linearization_point into
    // the vector linearized_stokes_initial_guess. We need to do the copy because
    // linearized_stokes_variables has a different
    // layout than current_linearization_point, which also contains all the
    // other solution variables.
    if (sim.assemble_newton_stokes_system == false)
      {
        linearized_stokes_initial_guess.block (block_vel) = sim.current_linearization_point.block (block_vel);
        linearized_stokes_initial_guess.block (block_p) = sim.current_linearization_point.block (block_p);

        sim.denormalize_pressure (sim.last_pressure_normalization_adjustment,
                                  linearized_stokes_initial_guess);
      }
    else
      {
        // The Newton solver solves for updates to variables, for which our best guess is zero when
        // the it isn't the first nonlinear iteration. When it is the first nonlinear iteration, we
        // have to assemble the full (non-defect correction) Picard, to get the boundary conditions
        // right in combination with being able to use the initial guess optimally. So we may never
        // end up here when it is the first nonlinear iteration.
        Assert(sim.nonlinear_iteration != 0,
               ExcMessage ("The Newton solver may not be active in the first nonlinear iteration"));

        linearized_stokes_initial_guess.block (block_vel) = 0;
        linearized_stokes_initial_guess.block (block_p) = 0;
      }

    sim.current_constraints.set_zero (linearized_stokes_initial_guess);
    linearized_stokes_initial_guess.block (block_p) /= sim.pressure_scaling;

    double solver_tolerance = 0;
    if (sim.assemble_newton_stokes_system == false)
      {
        // (ab)use the distributed solution vector to temporarily put a residual in
        // (we don't care about the residual vector -- all we care about is the
        // value (number) of the initial residual). The initial residual is returned
        // to the caller (for nonlinear computations). This value is computed before
        // the solve because we want to compute || A^{k+1} U^k - F^{k+1} ||, which is
        // the nonlinear residual. Because the place where the nonlinear residual is
        // checked against the nonlinear tolerance comes after the solve, the system
        // is solved one time too many in the case of a nonlinear Picard solver.

        // We must copy between Trilinos/dealii vector types
        dealii::LinearAlgebra::distributed::BlockVector<double> solution_copy(2);
        dealii::LinearAlgebra::distributed::BlockVector<double> initial_copy(2);
        dealii::LinearAlgebra::distributed::BlockVector<double> rhs_copy(2);

        stokes_matrix.initialize_dof_vector(solution_copy);
        stokes_matrix.initialize_dof_vector(initial_copy);
        stokes_matrix.initialize_dof_vector(rhs_copy);

        internal::ChangeVectorTypes::copy(solution_copy,distributed_stokes_solution);
        internal::ChangeVectorTypes::copy(initial_copy,linearized_stokes_initial_guess);
        internal::ChangeVectorTypes::copy(rhs_copy,distributed_stokes_rhs);

        // Compute residual l2_norm
        stokes_matrix.vmult(solution_copy,initial_copy);
        solution_copy.sadd(-1,1,rhs_copy);
        initial_nonlinear_residual = solution_copy.l2_norm();

        // Note: the residual is computed with a zero velocity, effectively computing
        // || B^T p - g ||, which we are going to use for our solver tolerance.
        // We do not use the current velocity for the initial residual because
        // this would not decrease the number of iterations if we had a better
        // initial guess (say using a smaller timestep). But we need to use
        // the pressure instead of only using the norm of the rhs, because we
        // are only interested in the part of the rhs not balanced by the static
        // pressure (the current pressure is a good approximation for the static
        // pressure).
        initial_copy.block(0) = 0.;
        stokes_matrix.vmult(solution_copy,initial_copy);
        solution_copy.block(0).sadd(-1,1,rhs_copy.block(0));

        const double residual_u = solution_copy.block(0).l2_norm();

        const double residual_p = rhs_copy.block(1).l2_norm();

        solver_tolerance = sim.parameters.linear_stokes_solver_tolerance *
                           std::sqrt(residual_u*residual_u+residual_p*residual_p);
      }
    else
      {
        // if we are solving for the Newton update, then the initial guess of the solution
        // vector is the zero vector, and the starting (nonlinear) residual is simply
        // the norm of the (Newton) right hand side vector
        const double residual_u = distributed_stokes_rhs.block(0).l2_norm();
        const double residual_p = distributed_stokes_rhs.block(1).l2_norm();
        solver_tolerance = sim.parameters.linear_stokes_solver_tolerance *
                           std::sqrt(residual_u*residual_u+residual_p*residual_p);

        // as described in the documentation of the function, the initial
        // nonlinear residual for the Newton method is computed by just
        // taking the norm of the right hand side
        initial_nonlinear_residual = std::sqrt(residual_u*residual_u+residual_p*residual_p);
      }

    // Now overwrite the solution vector again with the current best guess
    // to solve the linear system
    distributed_stokes_solution = linearized_stokes_initial_guess;

    // Again, copy solution and rhs vectors to solve with matrix-free operators
    dealii::LinearAlgebra::distributed::BlockVector<double> solution_copy(2);
    dealii::LinearAlgebra::distributed::BlockVector<double> rhs_copy(2);

    stokes_matrix.initialize_dof_vector(solution_copy);
    stokes_matrix.initialize_dof_vector(rhs_copy);

    internal::ChangeVectorTypes::copy(solution_copy,distributed_stokes_solution);
    internal::ChangeVectorTypes::copy(rhs_copy,distributed_stokes_rhs);

    // create Solver controls for the cheap and expensive solver phase
    SolverControl solver_control_cheap (sim.parameters.n_cheap_stokes_solver_steps,
                                        solver_tolerance, true);
    SolverControl solver_control_expensive (sim.parameters.n_expensive_stokes_solver_steps,
                                            solver_tolerance);

    solver_control_cheap.enable_history_data();
    solver_control_expensive.enable_history_data();

    // create a cheap preconditioner that consists of only a single V-cycle
    const internal::BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType, GMGPreconditioner, GMGPreconditioner>
    preconditioner_cheap (stokes_matrix, A_block_matrix, Schur_complement_block_matrix,
                          prec_A, prec_Schur,
                          /*do_solve_A*/false,
                          /*do_solve_Schur*/false,
                          sim.stokes_A_block_is_symmetric(),
                          sim.parameters.linear_solver_A_block_tolerance,
                          sim.parameters.linear_solver_S_block_tolerance);

    // create an expensive preconditioner that solves for the A block with CG
    const internal::BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType, GMGPreconditioner, GMGPreconditioner>
    preconditioner_expensive (stokes_matrix, A_block_matrix, Schur_complement_block_matrix,
                              prec_A, prec_Schur,
                              /*do_solve_A*/true,
                              /*do_solve_Schur*/true,
                              sim.stokes_A_block_is_symmetric(),
                              sim.parameters.linear_solver_A_block_tolerance,
                              sim.parameters.linear_solver_S_block_tolerance);

    PrimitiveVectorMemory<dealii::LinearAlgebra::distributed::BlockVector<double>> mem;

    // Time vmult of different matrix-free operators, solver IDR with the cheap preconditioner, and
    // solver GMRES with the cheap preconditioner. Each timing is repeated 10 times, and the
    // function may be called a couple of times within each timing, depending on the argument repeats.
    if (do_timings)
      {
        const int n_timings = 10;
        Timer timer(sim.mpi_communicator);

        auto time_this = [&](const char *name, int repeats, const std::function<void()> &body, const std::function<void()> &prepare)
        {
          sim.pcout << "Timing " << name << ' ' << n_timings << " time(s) and repeat "
                    << repeats << " time(s) within each timing:" << std::endl;

          body(); // warm up

          double average_time = 0.;

          for (int i=0; i<n_timings; ++i)
            {
              prepare();
              sim.pcout << "\t... " << std::flush;
              timer.restart();

              for (int r=0; r<repeats; ++r)
                body();

              timer.stop();
              double time = timer.wall_time();
              const double average_time_per_timing = time/repeats;
              sim.pcout << average_time_per_timing << std::endl;
              average_time += average_time_per_timing;
            }

          sim.pcout << "\taverage wall time of all: "<< average_time/n_timings << " seconds" << std::endl;

        };

        // stokes vmult
        {
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_dst = solution_copy;
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_src = rhs_copy;
          time_this("stokes_vmult", 10,
                    [&] ()
          {
            stokes_matrix.vmult(tmp_dst, tmp_src);
          },
          [&] ()
          {
            tmp_src = tmp_dst;
          }
                   );
        }

        // stokes preconditioner
        {
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_dst = solution_copy;
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_src = rhs_copy;
          time_this("stokes_preconditioner", 1,
                    [&] ()
          {
            preconditioner_cheap.vmult(tmp_dst, tmp_src);
          },
          [&] ()
          {
            tmp_src = tmp_dst;
          }
                   );
        }
        // A preconditioner
        {
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_dst = solution_copy;
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_src = rhs_copy;
          time_this("A_preconditioner", 1,
                    [&] ()
          {
            prec_A.vmult(tmp_dst.block(0), tmp_src.block(0));
          },
          [&] ()
          {
            tmp_src = tmp_dst;
          }
                   );
        }
        // S preconditioner
        {
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_dst = solution_copy;
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_src = rhs_copy;
          time_this("S_preconditioner", 5,
                    [&] ()
          {
            prec_Schur.vmult(tmp_dst.block(1), tmp_src.block(1));
          },
          [&] ()
          {
            tmp_src = tmp_dst;
          }
                   );
        }
        // Solve
        {
          // hard-code the number of iterations here to always do cheap iterations
          SolverControl solver_control_cheap (1000, solver_tolerance, true);

          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_dst = solution_copy;
          dealii::LinearAlgebra::distributed::BlockVector<double> tmp_src = rhs_copy;
          time_this("Stokes_solve_cheap_idr", 1,
                    [&]
          {
            SolverIDR<dealii::LinearAlgebra::distributed::BlockVector<double>>
            solver(solver_control_cheap, mem,
            SolverIDR<dealii::LinearAlgebra::distributed::BlockVector<double>>::
            AdditionalData(sim.parameters.idr_s_parameter));

            solver.solve (stokes_matrix,
            tmp_dst,
            tmp_src,
            preconditioner_cheap);
          },
          [&] ()
          {
            tmp_dst = solution_copy;
          }
                   );

          time_this("Stokes_solve_cheap_gmres", 1,
                    [&]
          {
            SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double>>
            solver(solver_control_cheap, mem,
            SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double>>::
            AdditionalData(sim.parameters.stokes_gmres_restart_length+2,
            true));

            solver.solve (stokes_matrix,
            tmp_dst,
            tmp_src,
            preconditioner_cheap);
          },
          [&] ()
          {
            tmp_dst = solution_copy;
          }
                   );
        }
      }

    // step 1a: try if the simple and fast solver
    // succeeds in n_cheap_stokes_solver_steps steps or less.
    try
      {
        // if this cheaper solver is not desired, then simply short-cut
        // the attempt at solving with the cheaper preconditioner
        if (sim.parameters.n_cheap_stokes_solver_steps == 0)
          throw SolverControl::NoConvergence(0,0);

        // Unlike with the expensive preconditioner which uses CG solves on both the
        // velocity and pressure space, the cheap preconditioner only contains matrix-vector
        // products and GMG v-cycle where the smoothers, transfer operators, and coarse
        // solvers are all defined to be linear operators which do not change from iteration
        // to iteration. Therefore we can use non-flexible Krylov methods like GMRES or IDR(s),
        // instead of requiring FGMRES, greatly lowing the memory requirement of the solver.
        if (sim.parameters.stokes_krylov_type == Parameters<dim>::StokesKrylovType::gmres)
          {
            SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double>>
            solver(solver_control_cheap, mem,
                   SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double>>::
                   AdditionalData(sim.parameters.stokes_gmres_restart_length+2,
                                  true));

            solver.solve (stokes_matrix,
                          solution_copy,
                          rhs_copy,
                          preconditioner_cheap);
          }
        else if (sim.parameters.stokes_krylov_type == Parameters<dim>::StokesKrylovType::idr_s)
          {
            SolverIDR<dealii::LinearAlgebra::distributed::BlockVector<double>>
            solver(solver_control_cheap, mem,
                   SolverIDR<dealii::LinearAlgebra::distributed::BlockVector<double>>::
                   AdditionalData(sim.parameters.idr_s_parameter));

            solver.solve (stokes_matrix,
                          solution_copy,
                          rhs_copy,
                          preconditioner_cheap);
          }
        else
          Assert(false,ExcNotImplemented());

        // Success. Print all iterations to screen (0 expensive iterations).
        sim.pcout << (solver_control_cheap.last_step() != numbers::invalid_unsigned_int ?
                      solver_control_cheap.last_step():
                      0)
                  << "+0"
                  << " iterations." << std::endl;

        final_linear_residual = solver_control_cheap.last_value();
      }
    // step 1b: take the stronger solver in case
    // the simple solver failed and attempt solving
    // it in n_expensive_stokes_solver_steps steps or less.
    catch (const SolverControl::NoConvergence &exc)
      {
        // The cheap solver failed or never ran.
        // Print the number of cheap iterations to screen to indicate we
        // try the expensive solver next.
        sim.pcout << (solver_control_cheap.last_step() != numbers::invalid_unsigned_int ?
                      solver_control_cheap.last_step():
                      0) << '+' << std::flush;

        // use the value defined by the user
        // OR
        // at least a restart length of 100 for melt models
        const unsigned int number_of_temporary_vectors = (sim.parameters.include_melt_transport == false ?
                                                          sim.parameters.stokes_gmres_restart_length :
                                                          std::max(sim.parameters.stokes_gmres_restart_length, 100U));

        SolverFGMRES<dealii::LinearAlgebra::distributed::BlockVector<double>>
        solver(solver_control_expensive, mem,
               SolverFGMRES<dealii::LinearAlgebra::distributed::BlockVector<double>>::
               AdditionalData(number_of_temporary_vectors));

        try
          {
            // if no expensive steps allowed, we have failed
            if (sim.parameters.n_expensive_stokes_solver_steps == 0)
              {
                sim.pcout << "0 iterations." << std::endl;
                throw exc;
              }

            solver.solve(stokes_matrix,
                         solution_copy,
                         rhs_copy,
                         preconditioner_expensive);

            // Success. Print expensive iterations to screen.
            sim.pcout << solver_control_expensive.last_step()
                      << " iterations." << std::endl;

            final_linear_residual = solver_control_expensive.last_value();
          }
        // if the solver fails trigger the post stokes solver signal and throw an exception
        catch (const std::exception &exc)
          {
            sim.signals.post_stokes_solver(sim,
                                           preconditioner_cheap.n_iterations_Schur_complement() + preconditioner_expensive.n_iterations_Schur_complement(),
                                           preconditioner_cheap.n_iterations_A_block() + preconditioner_expensive.n_iterations_A_block(),
                                           solver_control_cheap,
                                           solver_control_expensive);

            std::vector<SolverControl> solver_controls;
            if (sim.parameters.n_cheap_stokes_solver_steps > 0)
              solver_controls.push_back(solver_control_cheap);
            if (sim.parameters.n_expensive_stokes_solver_steps > 0)
              solver_controls.push_back(solver_control_expensive);

            Utilities::throw_linear_solver_failure_exception("iterative Stokes solver",
                                                             "StokesMatrixFreeHandlerImplementation::solve",
                                                             solver_controls,
                                                             exc,
                                                             sim.mpi_communicator,
                                                             sim.parameters.output_directory+"solver_history.txt");
          }
      }

    //signal successful solver
    sim.signals.post_stokes_solver(sim,
                                   preconditioner_cheap.n_iterations_Schur_complement() + preconditioner_expensive.n_iterations_Schur_complement(),
                                   preconditioner_cheap.n_iterations_A_block() + preconditioner_expensive.n_iterations_A_block(),
                                   solver_control_cheap,
                                   solver_control_expensive);

    // distribute hanging node and other constraints
    solution_copy.update_ghost_values();
    internal::ChangeVectorTypes::copy(distributed_stokes_solution,solution_copy);

#if DEAL_II_VERSION_GTE(9,6,0)
    IndexSet stokes_dofs (sim.dof_handler.n_dofs());
    stokes_dofs.add_range (0, distributed_stokes_solution.size());
    const AffineConstraints<double> current_stokes_constraints
      = sim.current_constraints.get_view (stokes_dofs);
    current_stokes_constraints.distribute(distributed_stokes_solution);
#else
    sim.current_constraints.distribute(distributed_stokes_solution);
#endif

    // now rescale the pressure back to real physical units
    distributed_stokes_solution.block(block_p) *= sim.pressure_scaling;

    // then copy back the solution from the temporary (non-ghosted) vector
    // into the ghosted one with all solution components
    solution_vector.block(block_vel) = distributed_stokes_solution.block(block_vel);
    solution_vector.block(block_p) = distributed_stokes_solution.block(block_p);

    if (print_details)
      {
        sim.pcout << "    Schur complement preconditioner: " << preconditioner_cheap.n_iterations_Schur_complement()
                  << '+'
                  << preconditioner_expensive.n_iterations_Schur_complement()
                  << " iterations." << std::endl;
        sim.pcout << "    A block preconditioner: " << preconditioner_cheap.n_iterations_A_block()
                  << '+'
                  << preconditioner_expensive.n_iterations_A_block()
                  << " iterations." << std::endl;
      }

    // do some cleanup now that we have the solution
    sim.remove_nullspace(solution_vector, distributed_stokes_solution);
    if (sim.assemble_newton_stokes_system == false)
      sim.last_pressure_normalization_adjustment = sim.normalize_pressure(solution_vector);


    // convert melt pressures
    // TODO: We assert in the StokesMatrixFreeHandler constructor that we
    //       are not including melt transport.
    if (sim.parameters.include_melt_transport)
      sim.melt_handler->compute_melt_variables(sim.system_matrix,solution_vector,sim.system_rhs);


    return std::pair<double,double>(initial_nonlinear_residual,
                                    final_linear_residual);
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::setup_dofs()
  {
    // Periodic boundary conditions with hanging nodes on the boundary currently
    // cause the GMG not to converge. We catch this case early to provide the
    // user with a reasonable error message:
    {
      bool have_periodic_hanging_nodes = false;
      for (const auto &cell : sim.triangulation.active_cell_iterators())
        if (cell->is_locally_owned())
          for (const auto f : cell->face_indices())
            {
              if (cell->has_periodic_neighbor(f))
                {
                  const auto &neighbor = cell->periodic_neighbor(f);
                  // This way, we can only detect the case where the neighbor is coarser,
                  // but this is fine as the other owner covers that situation:
                  if (neighbor->level()<cell->level())
                    have_periodic_hanging_nodes = true;
                }
            }

      have_periodic_hanging_nodes = (dealii::Utilities::MPI::max(have_periodic_hanging_nodes ? 1 : 0, sim.triangulation.get_communicator())) == 1;
      AssertThrow(have_periodic_hanging_nodes==false, ExcNotImplemented());
    }

    // This vector will be refilled with the new MatrixFree objects below:
    matrix_free_objects.clear();

    // Velocity DoFHandler
    {
      dof_handler_v.clear();
      dof_handler_v.distribute_dofs(fe_v);

      DoFRenumbering::hierarchical(dof_handler_v);

#if DEAL_II_VERSION_GTE(9,7,0)
      const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler_v);
#else
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler_v, locally_relevant_dofs);
#endif

#if DEAL_II_VERSION_GTE(9,6,0)
      constraints_v.reinit(dof_handler_v.locally_owned_dofs(), locally_relevant_dofs);
#else
      constraints_v.reinit(locally_relevant_dofs);
#endif

      {
        const auto &pbs = sim.geometry_model->get_periodic_boundary_pairs();

        for (const auto &p: pbs)
          {
            DoFTools::make_periodicity_constraints(dof_handler_v,
                                                   p.first.first,  // first boundary id
                                                   p.first.second, // second boundary id
                                                   p.second,       // cartesian direction for translational symmetry
                                                   constraints_v);
          }
      }
      DoFTools::make_hanging_node_constraints (dof_handler_v, constraints_v);
      sim.compute_initial_velocity_boundary_constraints(constraints_v);
      sim.compute_current_velocity_boundary_constraints(constraints_v);


      VectorTools::compute_no_normal_flux_constraints (dof_handler_v,
                                                       /* first_vector_component= */
                                                       0,
                                                       sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators(),
                                                       constraints_v,
                                                       *sim.mapping);
      constraints_v.close ();
    }

    // Pressure DoFHandler
    {
      dof_handler_p.clear();
      dof_handler_p.distribute_dofs(fe_p);

      DoFRenumbering::hierarchical(dof_handler_p);

#if DEAL_II_VERSION_GTE(9,7,0)
      const IndexSet locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler_p);
#else
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs(dof_handler_p,
                                              locally_relevant_dofs);
#endif

      constraints_p.reinit(
#if DEAL_II_VERSION_GTE(9,6,0)
        dof_handler_p.locally_owned_dofs(),
#endif
        locally_relevant_dofs);
      {
        const auto &pbs = sim.geometry_model->get_periodic_boundary_pairs();

        for (const auto &p: pbs)
          {
            DoFTools::make_periodicity_constraints(dof_handler_p,
                                                   p.first.first,  // first boundary id
                                                   p.first.second, // second boundary id
                                                   p.second,       // cartesian direction for translational symmetry
                                                   constraints_p);
          }
      }
      DoFTools::make_hanging_node_constraints (dof_handler_p, constraints_p);
      constraints_p.close();
    }

    // Coefficient transfer objects
    {
      dof_handler_projection.clear();
      dof_handler_projection.distribute_dofs(fe_projection);

      DoFRenumbering::hierarchical(dof_handler_projection);
    }

    // Multigrid DoF setup
    {
      //Ablock GMG
      dof_handler_v.distribute_mg_dofs();

      mg_constrained_dofs_A_block.clear();
      mg_constrained_dofs_A_block.initialize(dof_handler_v);

      std::set<types::boundary_id> dirichlet_boundary = sim.boundary_velocity_manager.get_zero_boundary_velocity_indicators();
      for (const auto boundary_id: sim.boundary_velocity_manager.get_prescribed_boundary_velocity_indicators())
        {
          const ComponentMask component_mask = sim.boundary_velocity_manager.get_component_mask(boundary_id);

          if (component_mask != ComponentMask(sim.introspection.n_components, false))
            {
              ComponentMask velocity_mask(fe_v.n_components(), false);

              for (unsigned int i=0; i<dim; ++i)
                velocity_mask.set(i, component_mask[sim.introspection.component_indices.velocities[i]]);

              mg_constrained_dofs_A_block.make_zero_boundary_constraints(dof_handler_v, {boundary_id}, velocity_mask);
            }
          else
            {
              // no mask given: add at the end
              dirichlet_boundary.insert(boundary_id);
            }
        }

      // Unconditionally call this function, even if the set is empty. Otherwise, the data structure
      // for boundary indices will not be created (if mesh has no Dirichlet conditions).
      mg_constrained_dofs_A_block.make_zero_boundary_constraints(dof_handler_v, dirichlet_boundary);

      //Schur complement matrix GMG
      dof_handler_p.distribute_mg_dofs();

      mg_constrained_dofs_Schur_complement.clear();
      mg_constrained_dofs_Schur_complement.initialize(dof_handler_p);

      dof_handler_projection.distribute_mg_dofs();
    }

    // Setup the matrix-free operators
    std::shared_ptr<MatrixFree<dim,double>> matrix_free = std::make_shared<MatrixFree<dim,double>>();
    matrix_free_objects.push_back(matrix_free);

    // Matrixfree object
    {
      typename MatrixFree<dim,double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_gradients | update_JxW_values);

      if (sim.mesh_deformation
          && !sim.mesh_deformation->get_free_surface_boundary_indicators().empty())
        additional_data.mapping_update_flags_boundary_faces =
          (update_values  |
           update_quadrature_points |
           update_normal_vectors |
           update_JxW_values);

      std::vector<const DoFHandler<dim>*> stokes_dofs {&dof_handler_v, &dof_handler_p};
      std::vector<const AffineConstraints<double> *> stokes_constraints {&constraints_v, &constraints_p};

      matrix_free->reinit(*sim.mapping, stokes_dofs, stokes_constraints,
                          QGauss<1>(sim.parameters.stokes_velocity_degree+1), additional_data);
    }

    // Stokes matrix
    {
      stokes_matrix.clear();
      stokes_matrix.initialize(matrix_free);
    }

    // ABlock matrix
    {
      A_block_matrix.clear();
      std::vector<unsigned int> selected = {0}; // select velocity DoFHandler
      A_block_matrix.initialize(matrix_free, selected);
    }

    // Schur complement block matrix
    {
      Schur_complement_block_matrix.clear();
      std::vector< unsigned int > selected = {1}; // select pressure DoFHandler
      Schur_complement_block_matrix.initialize(matrix_free, selected , selected);
    }

    // GMG matrices
    {
      const unsigned int n_levels = sim.triangulation.n_global_levels();

      mg_matrices_Schur_complement.clear_elements();
      mg_matrices_Schur_complement.resize(0, n_levels-1);
      mg_matrices_A_block.clear_elements();
      mg_matrices_A_block.resize(0, n_levels-1);

      for (unsigned int level=0; level<n_levels; ++level)
        {
          AffineConstraints<double> level_constraints_v;
          AffineConstraints<double> level_constraints_p;
          const Mapping<dim> &mapping =
            (sim.mesh_deformation) ? sim.mesh_deformation->get_level_mapping(level) : *sim.mapping;

          {
#if DEAL_II_VERSION_GTE(9,7,0)
            const IndexSet relevant_dofs = DoFTools::extract_locally_relevant_level_dofs(dof_handler_v, level);
#else
            IndexSet relevant_dofs;
            DoFTools::extract_locally_relevant_level_dofs(dof_handler_v, level, relevant_dofs);
#endif

#if DEAL_II_VERSION_GTE(9,6,0)
            level_constraints_v.reinit(dof_handler_v.locally_owned_mg_dofs(level), relevant_dofs);
            for (const auto index : mg_constrained_dofs_A_block.get_boundary_indices(level))
              level_constraints_v.constrain_dof_to_zero(index);
#else
            level_constraints_v.reinit(relevant_dofs);
            level_constraints_v.add_lines(mg_constrained_dofs_A_block.get_boundary_indices(level));
#endif
            level_constraints_v.close();

            std::set<types::boundary_id> no_flux_boundary
              = sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
            if (!no_flux_boundary.empty())
              {
                AffineConstraints<double> user_level_constraints;
#if DEAL_II_VERSION_GTE(9,6,0)
                user_level_constraints.reinit(dof_handler_v.locally_owned_mg_dofs(level), relevant_dofs);
#else
                user_level_constraints.reinit(relevant_dofs);
#endif
                const IndexSet &refinement_edge_indices =
                  mg_constrained_dofs_A_block.get_refinement_edge_indices(level);
                dealii::VectorTools::compute_no_normal_flux_constraints_on_level(
                  dof_handler_v,
                  0,
                  no_flux_boundary,
                  user_level_constraints,
                  mapping,
                  refinement_edge_indices,
                  level);

                user_level_constraints.close();
                mg_constrained_dofs_A_block.add_user_constraints(level,user_level_constraints);

                // let Dirichlet values win over no normal flux:
                level_constraints_v.merge(user_level_constraints, AffineConstraints<double>::left_object_wins);
                level_constraints_v.close();
              }
          }
          {
#if DEAL_II_VERSION_GTE(9,7,0)
            const IndexSet relevant_dofs = DoFTools::extract_locally_relevant_level_dofs(dof_handler_p, level);
#else
            IndexSet relevant_dofs;
            DoFTools::extract_locally_relevant_level_dofs(dof_handler_p, level, relevant_dofs);
#endif

#if DEAL_II_VERSION_GTE(9,6,0)
            level_constraints_p.reinit(dof_handler_p.locally_owned_mg_dofs(level), relevant_dofs);
#else
            level_constraints_p.reinit(relevant_dofs);
#endif

            level_constraints_p.close();
          }

          std::shared_ptr<MatrixFree<dim,GMGNumberType>> matrix_free_level = std::make_shared<MatrixFree<dim,GMGNumberType>>();
          matrix_free_objects.push_back(matrix_free_level);

          {
            typename MatrixFree<dim,GMGNumberType>::AdditionalData additional_data;
            additional_data.tasks_parallel_scheme = MatrixFree<dim,GMGNumberType>::AdditionalData::none;
            additional_data.mapping_update_flags = (update_gradients | update_JxW_values);
            additional_data.mg_level = level;

            std::vector<const DoFHandler<dim>*> stokes_dofs {&dof_handler_v, &dof_handler_p};
            std::vector<const AffineConstraints<double> *> stokes_constraints {&level_constraints_v,&level_constraints_p};

            matrix_free_level->reinit(mapping,
                                      stokes_dofs, stokes_constraints,
                                      QGauss<1>(sim.parameters.stokes_velocity_degree+1),
                                      additional_data);
          }
          {
            mg_matrices_A_block[level].clear();
            std::vector<unsigned int> selected = {0}; // select velocity DoFHandler
            mg_matrices_A_block[level].initialize(matrix_free_level, mg_constrained_dofs_A_block, level, selected);
          }
          {
            mg_matrices_Schur_complement[level].clear();
            std::vector<unsigned int> selected = {1}; // select pressure DoFHandler
            mg_matrices_Schur_complement[level].initialize(matrix_free_level, mg_constrained_dofs_Schur_complement, level, selected);
          }
        }
    }

    // Build MG transfer
    mg_transfer_A_block.clear();
    mg_transfer_A_block.initialize_constraints(mg_constrained_dofs_A_block);
    mg_transfer_A_block.build(dof_handler_v);

    mg_transfer_Schur_complement.clear();
    mg_transfer_Schur_complement.initialize_constraints(mg_constrained_dofs_Schur_complement);
    mg_transfer_Schur_complement.build(dof_handler_p);
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::build_preconditioner()
  {
    TimerOutput::Scope timer (this->sim.computing_timer, "Build Stokes preconditioner");

    for (unsigned int level=0; level < sim.triangulation.n_global_levels(); ++level)
      {
        mg_matrices_Schur_complement[level].compute_diagonal();
        mg_matrices_A_block[level].compute_diagonal();
      }
  }



  template <int dim, int velocity_degree>
  const DoFHandler<dim> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_dof_handler_v () const
  {
    return dof_handler_v;
  }



  template <int dim, int velocity_degree>
  const DoFHandler<dim> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_dof_handler_p () const
  {
    return dof_handler_p;
  }



  template <int dim, int velocity_degree>
  const DoFHandler<dim> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_dof_handler_projection () const
  {
    return dof_handler_projection;
  }



  template <int dim, int velocity_degree>
  const AffineConstraints<double> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_constraints_v() const
  {
    return constraints_v;
  }



  template <int dim, int velocity_degree>
  const AffineConstraints<double> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_constraints_p() const
  {
    return constraints_p;
  }



  template <int dim, int velocity_degree>
  const MGTransferMF<dim,GMGNumberType> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_mg_transfer_A() const
  {
    return mg_transfer_A_block;
  }



  template <int dim, int velocity_degree>
  const MGTransferMF<dim,GMGNumberType> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_mg_transfer_S() const
  {
    return mg_transfer_Schur_complement;
  }



  template <int dim, int velocity_degree>
  std::size_t
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>:: get_cell_data_memory_consumption() const
  {
    std::size_t total = active_cell_data.memory_consumption();

    for (unsigned int level=0; level<level_cell_data.max_level(); ++level)
      total += level_cell_data[level].memory_consumption();

    return total;
  }



// explicit instantiation of the functions we implement in this file
#define INSTANTIATE(dim) \
  template class StokesMatrixFreeHandler<dim>; \
  template class StokesMatrixFreeHandlerImplementation<dim,2>; \
  template class StokesMatrixFreeHandlerImplementation<dim,3>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
