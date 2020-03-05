/*
  Copyright (C) 2018 - 2019 by the authors of the ASPECT code.

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
#include <aspect/citation_info.h>
#include <aspect/melt.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/read_write_vector.templates.h>



namespace aspect
{
  namespace internal
  {

    /**
     * Here we define the function(s) to make no normal flux boundary constraints for
     * MG levels.
     */
    namespace TangentialBoundaryFunctions
    {
      template <int dim>
      void
      add_constraint(const std::array<types::global_dof_index,dim> &dof_indices,
                     const Tensor<1, dim> &constraining_vector,
                     ConstraintMatrix &constraints,
                     const double inhomogeneity = 0)
      {
        // This function is modified from an internal deal.II function in vector_tools.templates.h
        switch (dim)
          {
            case 2:
            {
              if (std::fabs(constraining_vector[0]) >
                  std::fabs(constraining_vector[1]) + 1e-10)
                {
                  if (!constraints.is_constrained(dof_indices[0]) &&
                      constraints.can_store_line(dof_indices[0]))
                    {
                      constraints.add_line(dof_indices[0]);

                      if (std::fabs(constraining_vector[1] /
                                    constraining_vector[0]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[0],
                                              dof_indices[1],
                                              -constraining_vector[1] /
                                              constraining_vector[0]);

                      if (std::fabs(inhomogeneity / constraining_vector[0]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.set_inhomogeneity(
                          dof_indices[0],
                          inhomogeneity / constraining_vector[0]);
                    }
                }
              else
                {
                  if (!constraints.is_constrained(dof_indices[1]) &&
                      constraints.can_store_line(dof_indices[1]))
                    {
                      constraints.add_line(dof_indices[1]);

                      if (std::fabs(constraining_vector[0] /
                                    constraining_vector[1]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[1],
                                              dof_indices[0],
                                              -constraining_vector[0] /
                                              constraining_vector[1]);

                      if (std::fabs(inhomogeneity / constraining_vector[1]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.set_inhomogeneity(
                          dof_indices[1],
                          inhomogeneity / constraining_vector[1]);
                    }
                }
              break;
            }

            case 3:
            {
              if ((std::fabs(constraining_vector[0]) >=
                   std::fabs(constraining_vector[1]) + 1e-10) &&
                  (std::fabs(constraining_vector[0]) >=
                   std::fabs(constraining_vector[2]) + 2e-10))
                {
                  if (!constraints.is_constrained(dof_indices[0]) &&
                      constraints.can_store_line(dof_indices[0]))
                    {
                      constraints.add_line(dof_indices[0]);

                      if (std::fabs(constraining_vector[1] /
                                    constraining_vector[0]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[0],
                                              dof_indices[1],
                                              -constraining_vector[1] /
                                              constraining_vector[0]);

                      if (std::fabs(constraining_vector[2] /
                                    constraining_vector[0]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[0],
                                              dof_indices[2],
                                              -constraining_vector[2] /
                                              constraining_vector[0]);

                      if (std::fabs(inhomogeneity / constraining_vector[0]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.set_inhomogeneity(
                          dof_indices[0],
                          inhomogeneity / constraining_vector[0]);
                    }
                }
              else if ((std::fabs(constraining_vector[1]) + 1e-10 >=
                        std::fabs(constraining_vector[0])) &&
                       (std::fabs(constraining_vector[1]) >=
                        std::fabs(constraining_vector[2]) + 1e-10))
                {
                  if (!constraints.is_constrained(dof_indices[1]) &&
                      constraints.can_store_line(dof_indices[1]))
                    {
                      constraints.add_line(dof_indices[1]);

                      if (std::fabs(constraining_vector[0] /
                                    constraining_vector[1]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[1],
                                              dof_indices[0],
                                              -constraining_vector[0] /
                                              constraining_vector[1]);

                      if (std::fabs(constraining_vector[2] /
                                    constraining_vector[1]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[1],
                                              dof_indices[2],
                                              -constraining_vector[2] /
                                              constraining_vector[1]);

                      if (std::fabs(inhomogeneity / constraining_vector[1]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.set_inhomogeneity(
                          dof_indices[1],
                          inhomogeneity / constraining_vector[1]);
                    }
                }
              else
                {
                  if (!constraints.is_constrained(dof_indices[2]) &&
                      constraints.can_store_line(dof_indices[2]))
                    {
                      constraints.add_line(dof_indices[2]);

                      if (std::fabs(constraining_vector[0] /
                                    constraining_vector[2]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[2],
                                              dof_indices[0],
                                              -constraining_vector[0] /
                                              constraining_vector[2]);

                      if (std::fabs(constraining_vector[1] /
                                    constraining_vector[2]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.add_entry(dof_indices[2],
                                              dof_indices[1],
                                              -constraining_vector[1] /
                                              constraining_vector[2]);

                      if (std::fabs(inhomogeneity / constraining_vector[2]) >
                          std::numeric_limits<double>::epsilon())
                        constraints.set_inhomogeneity(
                          dof_indices[2],
                          inhomogeneity / constraining_vector[2]);
                    }
                }

              break;
            }

            default:
              Assert(false, ExcNotImplemented());
          }
      }


      template <int dim, int spacedim>
      void compute_no_normal_flux_constraints_shell(const DoFHandler<dim,spacedim> &dof_handler,
                                                    const MGConstrainedDoFs        &mg_constrained_dofs,
                                                    const Mapping<dim> &mapping,
                                                    const unsigned int level,
                                                    const unsigned int first_vector_component,
                                                    const std::set<types::boundary_id> &boundary_ids,
                                                    ConstraintMatrix &constraints)
      {
        // TODO: This is a simplification of compute_no_normal_flux_constraints() from deal.II.
        // The differences are:
        // - It works on a specific level so we can ignore hanging nodes
        // - We use the normal vector given by the manifold (instead of averaging surface vectors)
        //
        // This should go into deal.II at some point, but it is too specific right now.

        const IndexSet &refinement_edge_indices = mg_constrained_dofs.get_refinement_edge_indices(level);

        const auto &fe = dof_handler.get_fe();
        const std::vector<Point<dim - 1>> &unit_support_points = fe.get_unit_face_support_points();
        const Quadrature<dim - 1> quadrature(unit_support_points);
        const unsigned int dofs_per_face = fe.dofs_per_face;
        std::vector<types::global_dof_index> face_dofs(dofs_per_face);


        FEFaceValues<dim, spacedim> fe_face_values(mapping,
                                                   fe,
                                                   quadrature,
                                                   update_quadrature_points |
                                                   update_normal_vectors);

        std::set<types::boundary_id>::iterator b_id;
        for (const auto &cell : dof_handler.cell_iterators_on_level(level))
          if (cell->level_subdomain_id() != numbers::artificial_subdomain_id
              &&
              cell->level_subdomain_id() != numbers::invalid_subdomain_id)
            for (unsigned int face_no = 0;
                 face_no < GeometryInfo<dim>::faces_per_cell;
                 ++face_no)
              if ((b_id = boundary_ids.find(cell->face(face_no)->boundary_id())) !=
                  boundary_ids.end())
                {
                  typename DoFHandler<dim, spacedim>::level_face_iterator face = cell->face(face_no);
                  face->get_mg_dof_indices(level, face_dofs);
                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int i = 0; i < face_dofs.size(); ++i)
                    if (fe.face_system_to_component_index(i).first ==
                        first_vector_component)
                      // Refinement edge indices are going to be constrained to 0 during a
                      // multigrid cycle and do not need no-normal-flux constraints, so skip them:
                      if (!refinement_edge_indices.is_element(face_dofs[i]))
                        {
                          const Point<dim> position = fe_face_values.quadrature_point(i);
                          std::array<types::global_dof_index,dim> dof_indices;
                          dof_indices[0] = face_dofs[i];
                          for (unsigned int k = 0; k < dofs_per_face; ++k)
                            if ((k != i) &&
                                (quadrature.point(k) == quadrature.point(i)) &&
                                (fe.face_system_to_component_index(k).first >=
                                 first_vector_component) &&
                                (fe.face_system_to_component_index(k).first <
                                 first_vector_component + dim))
                              dof_indices
                              [fe.face_system_to_component_index(k).first -
                               first_vector_component] = face_dofs[k];

                          Tensor<1, dim> normal_vector =
                            cell->face(face_no)->get_manifold().normal_vector(
                              cell->face(face_no), position);

                          // remove small entries:
                          for (unsigned int d = 0; d < dim; ++d)
                            if (std::fabs(normal_vector[d]) < 1e-13)
                              normal_vector[d] = 0;
                          normal_vector /= normal_vector.norm();

                          add_constraint<dim>(dof_indices, normal_vector, constraints, 0.0);
                        }
                }
      }

      template <int dim>
      void compute_no_normal_flux_constraints_box (const DoFHandler<dim>    &dof,
                                                   const types::boundary_id  bid,
                                                   const unsigned int first_vector_component,
                                                   MGConstrainedDoFs         &mg_constrained_dofs)
      {
        // For a given boundary id, find which vector component is on the boundary
        // and set a zero boundary constraint for those degrees of freedom.
        std::set<types::boundary_id> bid_set;
        bid_set.insert(bid);

        const unsigned int n_components = dof.get_fe_collection().n_components();
        Assert(first_vector_component + dim <= n_components,
               ExcIndexRange(first_vector_component, 0, n_components - dim + 1));

        ComponentMask comp_mask(n_components, false);


        typename Triangulation<dim>::face_iterator
        face = dof.get_triangulation().begin_face(),
        endf = dof.get_triangulation().end_face();
        for (; face != endf; ++face)
          if (face->boundary_id() == bid)
            for (unsigned int d = 0; d < dim; ++d)
              {
                Tensor<1, dim, double> unit_vec;
                unit_vec[d] = 1.0;

                Tensor<1, dim> normal_vec =
                  face->get_manifold().normal_vector(face, face->center());

                if (std::abs(std::abs(unit_vec * normal_vec) - 1.0) < 1e-10)
                  comp_mask.set(d + first_vector_component, true);
                else
                  Assert(
                    std::abs(unit_vec * normal_vec) < 1e-10,
                    ExcMessage(
                      "We can currently only support no normal flux conditions "
                      "for a specific boundary id if all faces are normal to the "
                      "x, y, or z axis."));
              }

        Assert(comp_mask.n_selected_components() == 1,
               ExcMessage(
                 "We can currently only support no normal flux conditions "
                 "for a specific boundary id if all faces are facing in the "
                 "same direction, i.e., a boundary normal to the x-axis must "
                 "have a different boundary id than a boundary normal to the "
                 "y- or z-axis and so on. If the mesh here was produced using "
                 "GridGenerator::..., setting colorize=true during mesh generation "
                 "and calling make_no_normal_flux_constraints() for each no normal "
                 "flux boundary will fulfill the condition."));

        mg_constrained_dofs.make_zero_boundary_constraints(dof, bid_set, comp_mask);
      }
    }

    /**
     * Matrix-free operators must use deal.II defined vectors, while the rest of the ASPECT
     * software is based on Trilinos vectors. Here we define functions which copy between the
     * vector types.
     */
    namespace ChangeVectorTypes
    {
      void import(TrilinosWrappers::MPI::Vector &out,
                  const dealii::LinearAlgebra::ReadWriteVector<double> &rwv,
                  const VectorOperation::values                 operation)
      {
        Assert(out.size() == rwv.size(),
               ExcMessage("Both vectors need to have the same size for import() to work!"));

        Assert(out.locally_owned_elements() == rwv.get_stored_elements(),
               ExcNotImplemented());

        if (operation == VectorOperation::insert)
          {
            for (const auto idx : out.locally_owned_elements())
              out[idx] = rwv[idx];
          }
        else if (operation == VectorOperation::add)
          {
            for (const auto idx : out.locally_owned_elements())
              out[idx] += rwv[idx];
          }
        else
          AssertThrow(false, ExcNotImplemented());

        out.compress(operation);
      }


      void copy(TrilinosWrappers::MPI::Vector &out,
                const dealii::LinearAlgebra::distributed::Vector<double> &in)
      {
        dealii::LinearAlgebra::ReadWriteVector<double> rwv(out.locally_owned_elements());
        rwv.import(in, VectorOperation::insert);
        //This import function doesn't exist until after dealii 9.0
        //Implemented above
        import(out, rwv,VectorOperation::insert);
      }

      void copy(dealii::LinearAlgebra::distributed::Vector<double> &out,
                const TrilinosWrappers::MPI::Vector &in)
      {
        dealii::LinearAlgebra::ReadWriteVector<double> rwv;
        rwv.reinit(in);
        out.import(rwv, VectorOperation::insert);
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
        mutable unsigned int                                            n_iterations_A_;
        mutable unsigned int                                            n_iterations_Schur_complement_;
        const double                                                    A_block_tolerance;
        const double                                                    Schur_complement_tolerance;
        mutable dealii::LinearAlgebra::distributed::BlockVector<double> utmp;
        mutable dealii::LinearAlgebra::distributed::BlockVector<double> ptmp;
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
        {
          utmp.reinit(src);
          ptmp.reinit(src);
        }

      // either solve with the Schur complement matrix (if do_solve_Schur_complement==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_Schur_complement)
        {
          // first solve with the bottom left block, which we have built
          // as a mass matrix with the inverse of the viscosity
          SolverControl solver_control(100, src.block(1).l2_norm() * Schur_complement_tolerance,true);

          SolverCG<dealii::LinearAlgebra::distributed::Vector<double> > solver(solver_control);
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
                  dst.block(1) = 0.0;
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
                  if (Utilities::MPI::this_mpi_process(src.block(0).get_mpi_communicator()) == 0)
                    AssertThrow (false,
                                 ExcMessage (std::string("The iterative (bottom right) solver in BlockSchurGMGPreconditioner::vmult "
                                                         "did not converge to a tolerance of "
                                                         + Utilities::to_string(solver_control.tolerance()) +
                                                         ". It reported the following error:\n\n")
                                             +
                                             exc.what()))
                    else
                      throw QuietException();
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
        ptmp = dst;
        ptmp.block(0) = 0.0;
        stokes_matrix.vmult(utmp, ptmp); // B^T
        utmp.block(0) *= -1.0;
        utmp.block(0) += src.block(0);
      }

      // now either solve with the top left block (if do_solve_A==true)
      // or just apply one preconditioner sweep (for the first few
      // iterations of our two-stage outer GMRES iteration)
      if (do_solve_A == true)
        {
          SolverControl solver_control(1000, utmp.block(0).l2_norm() * A_block_tolerance);
          SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> solver(solver_control);
          try
            {
              dst.block(0) = 0.0;
              solver.solve(A_block, dst.block(0), utmp.block(0),
                           A_block_preconditioner);
              n_iterations_A_ += solver_control.last_step();
            }
          // if the solver fails, report the error from processor 0 with some additional
          // information about its location, and throw a quiet exception on all other
          // processors
          catch (const std::exception &exc)
            {
              if (Utilities::MPI::this_mpi_process(src.block(0).get_mpi_communicator()) == 0)
                AssertThrow (false,
                             ExcMessage (std::string("The iterative (top left) solver in BlockSchurGMGPreconditioner::vmult "
                                                     "did not converge to a tolerance of "
                                                     + Utilities::to_string(solver_control.tolerance()) +
                                                     ". It reported the following error:\n\n")
                                         +
                                         exc.what()))
                else
                  throw QuietException();
            }

        }
      else
        {
          A_block_preconditioner.vmult (dst.block(0), utmp.block(0));
          n_iterations_A_ += 1;
        }
    }
  }

  /**
   * Implementation of the matrix-free operators.
   *
   * Stokes operator
   */
  template <int dim, int degree_v, typename number>
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::StokesOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number> >()
  {}

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::clear ()
  {
    viscosity_x_2.reinit(TableIndices<1>(0));
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::BlockVector<number> >::clear();
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::
  fill_cell_data (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                  const double pressure_scaling,
                  const Triangulation<dim> &tria,
                  const DoFHandler<dim> &dof_handler_for_projection,
                  const bool is_compressible)
  {
    const unsigned int n_cells = this->data->n_macro_cells();
    viscosity_x_2.reinit(TableIndices<1>(n_cells));

    std::vector<types::global_dof_index> local_dof_indices(dof_handler_for_projection.get_fe().dofs_per_cell);
    for (unsigned int cell=0; cell<n_cells; ++cell)
      for (unsigned int i=0; i<this->get_matrix_free()->n_components_filled(cell); ++i)
        {
          typename DoFHandler<dim>::active_cell_iterator FEQ_cell = this->get_matrix_free()->get_cell_iterator(cell,i);
          typename DoFHandler<dim>::active_cell_iterator DG_cell(&tria,
                                                                 FEQ_cell->level(),
                                                                 FEQ_cell->index(),
                                                                 &dof_handler_for_projection);
          DG_cell->get_active_or_mg_dof_indices(local_dof_indices);

          //TODO: projection with higher degree
          Assert(local_dof_indices.size() == 1, ExcNotImplemented());
          viscosity_x_2(cell)[i] = 2.0*viscosity_values(local_dof_indices[0]);
        }

    this->pressure_scaling = pressure_scaling;
    this->is_compressible = is_compressible;
  }

  template <int dim, int degree_v, typename number>
  const Table<1, VectorizedArray<number> > &
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::get_viscosity_x_2_table()
  {
    return viscosity_x_2;
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>
  ::compute_diagonal ()
  {
    // There is currently no need in the code for the diagonal of the entire stokes
    // block. If needed, one could easily construct based on the diagonal of the A
    // block and append zeros to the end for the number of pressure DoFs.
    Assert(false, ExcNotImplemented());
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>
  ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                 dealii::LinearAlgebra::distributed::BlockVector<number>       &dst,
                 const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                 const std::pair<unsigned int, unsigned int>           &cell_range) const
  {
    FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (data, 0);
    FEEvaluation<dim,degree_v-1,  degree_v+1,1,  number> pressure (data, 1);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        const VectorizedArray<number> &cell_viscosity_x_2 = viscosity_x_2(cell);

        velocity.reinit (cell);
        velocity.read_dof_values (src.block(0));
        velocity.evaluate (false,true,false);
        pressure.reinit (cell);
        pressure.read_dof_values (src.block(1));
        pressure.evaluate (true,false,false);

        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
                                                          velocity.get_symmetric_gradient (q);
            VectorizedArray<number> pres = pressure.get_value(q);
            VectorizedArray<number> div = trace(sym_grad_u);
            pressure.submit_value(-1.0*pressure_scaling*div, q);

            sym_grad_u *= cell_viscosity_x_2;

            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= pressure_scaling*pres;

            if (is_compressible)
              for (unsigned int d=0; d<dim; ++d)
                sym_grad_u[d][d] -= cell_viscosity_x_2/3.0*div;

            velocity.submit_symmetric_gradient(sym_grad_u, q);
          }

        velocity.integrate (false,true);
        velocity.distribute_local_to_global (dst.block(0));
        pressure.integrate (true,false);
        pressure.distribute_local_to_global (dst.block(1));
      }
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>
  ::apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
               const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const
  {
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number> >::
    data->cell_loop(&StokesOperator::local_apply, this, dst, src);
  }

  /**
   * Mass matrix operator on pressure
   */
  template <int dim, int degree_p, typename number>
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::MassMatrixOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number> >()
  {}

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::clear ()
  {
    one_over_viscosity.reinit(TableIndices<1>(0));
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::clear();
  }

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::
  fill_cell_data (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                  const Triangulation<dim> &tria,
                  const DoFHandler<dim> &dof_handler_for_projection,
                  const bool is_mg_level_data,
                  const double pressure_scaling)
  {
    const unsigned int n_cells = this->data->n_macro_cells();
    one_over_viscosity.reinit(TableIndices<1>(n_cells));

    std::vector<types::global_dof_index> local_dof_indices(dof_handler_for_projection.get_fe().dofs_per_cell);
    for (unsigned int cell=0; cell<n_cells; ++cell)
      for (unsigned int i=0; i<this->get_matrix_free()->n_components_filled(cell); ++i)
        {
          if (is_mg_level_data)
            {
              typename DoFHandler<dim>::level_cell_iterator FEQ_cell = this->get_matrix_free()->get_cell_iterator(cell,i);
              typename DoFHandler<dim>::level_cell_iterator DG_cell(&tria,
                                                                    FEQ_cell->level(),
                                                                    FEQ_cell->index(),
                                                                    &dof_handler_for_projection);
              DG_cell->get_active_or_mg_dof_indices(local_dof_indices);
            }
          else
            {
              typename DoFHandler<dim>::active_cell_iterator FEQ_cell = this->get_matrix_free()->get_cell_iterator(cell,i);
              typename DoFHandler<dim>::active_cell_iterator DG_cell(&tria,
                                                                     FEQ_cell->level(),
                                                                     FEQ_cell->index(),
                                                                     &dof_handler_for_projection);
              DG_cell->get_active_or_mg_dof_indices(local_dof_indices);
            }

          //TODO: projection with higher degree
          Assert(local_dof_indices.size() == 1, ExcNotImplemented());
          one_over_viscosity(cell)[i] = 1.0/viscosity_values(local_dof_indices[0]);
        }

    this->pressure_scaling = pressure_scaling;
  }

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                 dealii::LinearAlgebra::distributed::Vector<number>       &dst,
                 const dealii::LinearAlgebra::distributed::Vector<number> &src,
                 const std::pair<unsigned int, unsigned int>           &cell_range) const
  {
    FEEvaluation<dim,degree_p,degree_p+2,1,number> pressure (data);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        const VectorizedArray<number> &cell_one_over_viscosity = one_over_viscosity(cell);

        pressure.reinit (cell);
        pressure.read_dof_values(src);
        pressure.evaluate (true, false);
        for (unsigned int q=0; q<pressure.n_q_points; ++q)
          pressure.submit_value(cell_one_over_viscosity*pressure_scaling*pressure_scaling*
                                pressure.get_value(q),q);
        pressure.integrate (true, false);
        pressure.distribute_local_to_global (dst);
      }
  }

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
               const dealii::LinearAlgebra::distributed::Vector<number> &src) const
  {
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::
    data->cell_loop(&MassMatrixOperator::local_apply, this, dst, src);
  }

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::compute_diagonal ()
  {
    this->inverse_diagonal_entries.
    reset(new DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number> >());
    this->diagonal_entries.
    reset(new DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number> >());

    dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    dealii::LinearAlgebra::distributed::Vector<number> &diagonal =
      this->diagonal_entries->get_vector();

    unsigned int dummy = 0;
    this->data->initialize_dof_vector(inverse_diagonal);
    this->data->initialize_dof_vector(diagonal);

    this->data->cell_loop (&MassMatrixOperator::local_compute_diagonal, this,
                           diagonal, dummy);

    this->set_constrained_entries_to_one(diagonal);
    inverse_diagonal = diagonal;
    const unsigned int local_size = inverse_diagonal.local_size();
    for (unsigned int i=0; i<local_size; ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i)
          =1./inverse_diagonal.local_element(i);
      }
  }

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>
  ::local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                            dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                            const unsigned int &,
                            const std::pair<unsigned int,unsigned int>       &cell_range) const
  {
    FEEvaluation<dim,degree_p,degree_p+2,1,number> pressure (data, 0);
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        const VectorizedArray<number> &cell_one_over_viscosity = one_over_viscosity(cell);

        pressure.reinit (cell);
        AlignedVector<VectorizedArray<number> > diagonal(pressure.dofs_per_cell);
        for (unsigned int i=0; i<pressure.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<pressure.dofs_per_cell; ++j)
              pressure.begin_dof_values()[j] = VectorizedArray<number>();
            pressure.begin_dof_values()[i] = make_vectorized_array<number> (1.);

            pressure.evaluate (true,false,false);
            for (unsigned int q=0; q<pressure.n_q_points; ++q)
              pressure.submit_value(cell_one_over_viscosity*pressure_scaling*pressure_scaling*
                                    pressure.get_value(q),q);
            pressure.integrate (true,false);

            diagonal[i] = pressure.begin_dof_values()[i];
          }

        for (unsigned int i=0; i<pressure.dofs_per_cell; ++i)
          pressure.begin_dof_values()[i] = diagonal[i];
        pressure.distribute_local_to_global (dst);
      }
  }

  /**
   * Velocity block operator
   */
  template <int dim, int degree_v, typename number>
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::ABlockOperator ()
    :
    MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number> >()
  {}

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::clear ()
  {
    viscosity_x_2.reinit(TableIndices<1>(0));
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::clear();
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::
  fill_cell_data (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                  const Triangulation<dim> &tria,
                  const DoFHandler<dim> &dof_handler_for_projection,
                  const bool is_mg_level_data,
                  const bool is_compressible)
  {
    const unsigned int n_cells = this->data->n_macro_cells();
    viscosity_x_2.reinit(TableIndices<1>(n_cells));

    std::vector<types::global_dof_index> local_dof_indices(dof_handler_for_projection.get_fe().dofs_per_cell);
    for (unsigned int cell=0; cell<n_cells; ++cell)
      for (unsigned int i=0; i<this->get_matrix_free()->n_components_filled(cell); ++i)
        {
          if (is_mg_level_data)
            {
              typename DoFHandler<dim>::level_cell_iterator FEQ_cell = this->get_matrix_free()->get_cell_iterator(cell,i);
              typename DoFHandler<dim>::level_cell_iterator DG_cell(&tria,
                                                                    FEQ_cell->level(),
                                                                    FEQ_cell->index(),
                                                                    &dof_handler_for_projection);
              DG_cell->get_active_or_mg_dof_indices(local_dof_indices);
            }
          else
            {
              typename DoFHandler<dim>::active_cell_iterator FEQ_cell = this->get_matrix_free()->get_cell_iterator(cell,i);
              typename DoFHandler<dim>::active_cell_iterator DG_cell(&tria,
                                                                     FEQ_cell->level(),
                                                                     FEQ_cell->index(),
                                                                     &dof_handler_for_projection);
              DG_cell->get_active_or_mg_dof_indices(local_dof_indices);
            }

          //TODO: projection with higher degree
          Assert(local_dof_indices.size() == 1, ExcNotImplemented());
          viscosity_x_2(cell)[i] = 2.0*viscosity_values(local_dof_indices[0]);
        }

    this->is_compressible = is_compressible;
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
        const VectorizedArray<number> &cell_viscosity_x_2 = viscosity_x_2(cell);

        velocity.reinit (cell);
        velocity.read_dof_values(src);
        velocity.evaluate (false, true, false);
        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
                                                          velocity.get_symmetric_gradient (q);
            sym_grad_u *= cell_viscosity_x_2;

            if (is_compressible)
              {
                VectorizedArray<number> div = trace(sym_grad_u);
                for (unsigned int d=0; d<dim; ++d)
                  sym_grad_u[d][d] -= 1.0/3.0*div;
              }
            velocity.submit_symmetric_gradient(sym_grad_u, q);
          }
        velocity.integrate (false, true);
        velocity.distribute_local_to_global (dst);
      }
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
               const dealii::LinearAlgebra::distributed::Vector<number> &src) const
  {
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::
    data->cell_loop(&ABlockOperator::local_apply, this, dst, src);
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::compute_diagonal ()
  {
    this->inverse_diagonal_entries.
    reset(new DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number> >());
    dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    unsigned int dummy = 0;
    this->data->cell_loop (&ABlockOperator::local_compute_diagonal, this,
                           inverse_diagonal, dummy);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i=0; i<inverse_diagonal.local_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1./inverse_diagonal.local_element(i);
      }
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                            dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                            const unsigned int &,
                            const std::pair<unsigned int,unsigned int>       &cell_range) const
  {
    FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (data, 0);
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        const VectorizedArray<number> &cell_viscosity_x_2 = viscosity_x_2(cell);

        velocity.reinit (cell);
        AlignedVector<VectorizedArray<number> > diagonal(velocity.dofs_per_cell);
        for (unsigned int i=0; i<velocity.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<velocity.dofs_per_cell; ++j)
              velocity.begin_dof_values()[j] = VectorizedArray<number>();
            velocity.begin_dof_values()[i] = make_vectorized_array<number> (1.);

            velocity.evaluate (false,true,false);
            for (unsigned int q=0; q<velocity.n_q_points; ++q)
              {
                SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
                                                              velocity.get_symmetric_gradient (q);

                sym_grad_u *= cell_viscosity_x_2;

                if (is_compressible)
                  {
                    VectorizedArray<number> div = trace(sym_grad_u);
                    for (unsigned int d=0; d<dim; ++d)
                      sym_grad_u[d][d] -= 1.0/3.0*div;
                  }

                velocity.submit_symmetric_gradient(sym_grad_u, q);
              }
            velocity.integrate (false,true);

            diagonal[i] = velocity.begin_dof_values()[i];
          }

        for (unsigned int i=0; i<velocity.dofs_per_cell; ++i)
          velocity.begin_dof_values()[i] = diagonal[i];
        velocity.distribute_local_to_global (dst);
      }
  }



  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>
  ::set_diagonal (const dealii::LinearAlgebra::distributed::Vector<number> &diag)
  {
    this->inverse_diagonal_entries.
    reset(new DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<number> >());
    dealii::LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);

    inverse_diagonal = diag;

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i=0; i<inverse_diagonal.local_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1./inverse_diagonal.local_element(i);
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

      stokes_fe (FE_Q<dim>(sim.parameters.stokes_velocity_degree),dim,
                 FE_Q<dim>(sim.parameters.stokes_velocity_degree-1),1),
      fe_v (FE_Q<dim>(sim.parameters.stokes_velocity_degree), dim),
      fe_p (FE_Q<dim>(sim.parameters.stokes_velocity_degree-1),1),
      fe_projection(FE_DGQ<dim>(0),1)
  {
    parse_parameters(prm);
    CitationInfo::add("mf");

    // This requires: porting the additional stabilization terms and using a
    // different mapping in the MatrixFree operators:
    AssertThrow(!sim.parameters.mesh_deformation_enabled, ExcNotImplemented());
    // Sorry, not any time soon:
    AssertThrow(!sim.parameters.include_melt_transport, ExcNotImplemented());
    // Not very difficult to do, but will require a different mass matrix
    // operator:
    AssertThrow(!sim.parameters.use_locally_conservative_discretization, ExcNotImplemented());


    // sanity check:
    Assert(sim.introspection.variable("velocity").block_index==0, ExcNotImplemented());
    Assert(sim.introspection.variable("pressure").block_index==1, ExcNotImplemented());

    // This is not terribly complicated, but we need to check that constraints
    // are set correctly, that the preconditioner converges, and requires
    // testing.
    AssertThrow(sim.geometry_model->get_periodic_boundary_pairs().size()==0, ExcNotImplemented());

    // We currently only support averaging that gives a constant value:
    using avg = MaterialModel::MaterialAveraging::AveragingOperation;
    AssertThrow((sim.parameters.material_averaging &
                 (avg::arithmetic_average | avg::harmonic_average | avg::geometric_average
                  | avg::pick_largest | avg::log_average))!=0
                , ExcMessage("The matrix-free Stokes solver currently only works if material model averaging is enabled"));

    // Currently cannot solve compressible flow with implicit reference density
    if (sim.material_model->is_compressible() == true)
      AssertThrow(sim.parameters.formulation_mass_conservation !=
                  Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile,
                  ExcNotImplemented());

    {
      const unsigned int n_vect_doubles =
        VectorizedArray<double>::n_array_elements;
      const unsigned int n_vect_bits = 8 * sizeof(double) * n_vect_doubles;

      sim.pcout << "Vectorization over " << n_vect_doubles
                << " doubles = " << n_vect_bits << " bits ("
                << dealii::Utilities::System::get_current_vectorization_level()
                << "), VECTORIZATION_LEVEL=" << DEAL_II_COMPILER_VECTORIZATION_LEVEL
                << std::endl;
    }
  }


  template <int dim, int velocity_degree>
  double StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_workload_imbalance ()
  {
    unsigned int n_proc = Utilities::MPI::n_mpi_processes(sim.triangulation.get_communicator());
    unsigned int n_global_levels = sim.triangulation.n_global_levels();

    unsigned long long int work_estimate = 0;
    unsigned long long int total_cells_in_hierarchy = 0;

    for (int lvl=n_global_levels-1; lvl>=0; --lvl)
      {
        unsigned long long int work_estimate_this_level;
        unsigned long long int total_cells_on_lvl;
        unsigned long long int n_owned_cells_on_lvl = 0;

        for (const auto &cell: sim.triangulation.cell_iterators_on_level(lvl))
          if (cell->is_locally_owned_on_level())
            n_owned_cells_on_lvl += 1;

        work_estimate_this_level = dealii::Utilities::MPI::max(n_owned_cells_on_lvl,sim.triangulation.get_communicator());

        work_estimate += work_estimate_this_level;

        total_cells_on_lvl = dealii::Utilities::MPI::sum(n_owned_cells_on_lvl,sim.triangulation.get_communicator());

        total_cells_in_hierarchy += total_cells_on_lvl;
      }
    double ideal_work = static_cast<double>(total_cells_in_hierarchy) / static_cast<double>(n_proc);
    double workload_imbalance_ratio = work_estimate / ideal_work;

    return workload_imbalance_ratio;
  }


  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::evaluate_material_model ()
  {
    {
      const QGauss<dim> quadrature_formula (sim.parameters.stokes_velocity_degree+1);

      FEValues<dim> fe_values (*sim.mapping,
                               sim.finite_element,
                               quadrature_formula,
                               update_values   |
                               update_gradients |
                               update_quadrature_points |
                               update_JxW_values);

      MaterialModel::MaterialModelInputs<dim> in(fe_values.n_quadrature_points, sim.introspection.n_compositional_fields);
      MaterialModel::MaterialModelOutputs<dim> out(fe_values.n_quadrature_points, sim.introspection.n_compositional_fields);

      std::vector<types::global_dof_index> local_dof_indices(fe_projection.dofs_per_cell);
      active_coef_dof_vec = 0.;

      // compute the integral quantities by quadrature
      for (const auto &cell: sim.dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            in.reinit(fe_values, cell, sim.introspection, sim.current_linearization_point);

            sim.material_model->fill_additional_material_model_inputs(in, sim.current_linearization_point, fe_values, sim.introspection);
            sim.material_model->evaluate(in, out);

            MaterialModel::MaterialAveraging::average (sim.parameters.material_averaging,
                                                       cell,
                                                       quadrature_formula,
                                                       *sim.mapping,
                                                       out);

            // we grab the first value, but all of them should be averaged to the same value:
            const double viscosity = out.viscosities[0];

            typename DoFHandler<dim>::active_cell_iterator dg_cell(&sim.triangulation,
                                                                   cell->level(),
                                                                   cell->index(),
                                                                   &dof_handler_projection);
            dg_cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < fe_projection.dofs_per_cell; ++i)
              active_coef_dof_vec[local_dof_indices[i]] = viscosity;
          }
      active_coef_dof_vec.compress(VectorOperation::insert);
    }

    const bool is_compressible = sim.material_model->is_compressible();

    stokes_matrix.fill_cell_data(active_coef_dof_vec,
                                 sim.pressure_scaling,
                                 sim.triangulation,
                                 dof_handler_projection,
                                 is_compressible);

    if (sim.parameters.n_expensive_stokes_solver_steps > 0)
      {
        A_block_matrix.fill_cell_data(active_coef_dof_vec,
                                      sim.triangulation,
                                      dof_handler_projection,
                                      /*is_mg_level_data*/false,
                                      is_compressible);

        Schur_complement_block_matrix.fill_cell_data(active_coef_dof_vec,
                                                     sim.triangulation,
                                                     dof_handler_projection,
                                                     /*is_mg_level_data*/false,
                                                     sim.pressure_scaling);
      }


    // Project to MG
    const unsigned int n_levels = sim.triangulation.n_global_levels();
    level_coef_dof_vec = 0.;
    level_coef_dof_vec.resize(0,n_levels-1);

    MGTransferMatrixFree<dim,double> transfer;
    transfer.build(dof_handler_projection);
    transfer.interpolate_to_mg(dof_handler_projection,
                               level_coef_dof_vec,
                               active_coef_dof_vec);

    for (unsigned int level=0; level<n_levels; ++level)
      {
        mg_matrices_A_block[level].fill_cell_data(level_coef_dof_vec[level],
                                                  sim.triangulation,
                                                  dof_handler_projection,
                                                  /*is_mg_level_data*/true,
                                                  is_compressible);

        mg_matrices_Schur_complement[level].fill_cell_data(level_coef_dof_vec[level],
                                                           sim.triangulation,
                                                           dof_handler_projection,
                                                           /*is_mg_level_data*/true,
                                                           sim.pressure_scaling);
      }
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::correct_stokes_rhs()
  {
    dealii::LinearAlgebra::distributed::BlockVector<double> rhs_correction(2);
    dealii::LinearAlgebra::distributed::BlockVector<double> u0(2);

    stokes_matrix.initialize_dof_vector(rhs_correction);
    stokes_matrix.initialize_dof_vector(u0);

    rhs_correction.collect_sizes();
    u0.collect_sizes();

    u0 = 0;
    rhs_correction = 0;
    sim.current_constraints.distribute(u0);
    u0.update_ghost_values();

    FEEvaluation<dim,velocity_degree,velocity_degree+1,dim,double>
    velocity (*stokes_matrix.get_matrix_free(), 0);
    FEEvaluation<dim,velocity_degree-1,velocity_degree+1,1,double>
    pressure (*stokes_matrix.get_matrix_free(), 1);

    for (unsigned int cell=0; cell<stokes_matrix.get_matrix_free()->n_macro_cells(); ++cell)
      {
        const VectorizedArray<double> &cell_viscosity_x_2 = stokes_matrix.get_viscosity_x_2_table()(cell);

        velocity.reinit (cell);
        velocity.read_dof_values_plain (u0.block(0));
        velocity.evaluate (false,true,false);
        pressure.reinit (cell);
        pressure.read_dof_values_plain (u0.block(1));
        pressure.evaluate (true,false,false);

        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            SymmetricTensor<2,dim,VectorizedArray<double>> sym_grad_u =
                                                          velocity.get_symmetric_gradient (q);
            VectorizedArray<double> pres = pressure.get_value(q);
            VectorizedArray<double> div = -trace(sym_grad_u);
            pressure.submit_value   (-1.0*sim.pressure_scaling*div, q);

            sym_grad_u *= cell_viscosity_x_2;

            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= sim.pressure_scaling*pres;

            velocity.submit_symmetric_gradient(-1.0*sym_grad_u, q);
          }

        velocity.integrate (false,true);
        velocity.distribute_local_to_global (rhs_correction.block(0));
        pressure.integrate (true,false);
        pressure.distribute_local_to_global (rhs_correction.block(1));
      }
    rhs_correction.compress(VectorOperation::add);

    LinearAlgebra::BlockVector stokes_rhs_correction (sim.introspection.index_sets.stokes_partitioning, sim.mpi_communicator);
    internal::ChangeVectorTypes::copy(stokes_rhs_correction,rhs_correction);
    sim.system_rhs.block(0) += stokes_rhs_correction.block(0);
    sim.system_rhs.block(1) += stokes_rhs_correction.block(1);
  }



  template <int dim, int velocity_degree>
  std::pair<double,double> StokesMatrixFreeHandlerImplementation<dim,velocity_degree>::solve()
  {
    double initial_nonlinear_residual = numbers::signaling_nan<double>();
    double final_linear_residual      = numbers::signaling_nan<double>();

    // Below we define all the objects needed to build the GMG preconditioner:
    using vector_t = dealii::LinearAlgebra::distributed::Vector<double>;

    // ABlock GMG Smoother: Chebyshev, degree 4
    typedef PreconditionChebyshev<ABlockMatrixType,vector_t> ASmootherType;
    mg::SmootherRelaxation<ASmootherType, vector_t>
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
              smoother_data_A[0].degree = numbers::invalid_unsigned_int;
              smoother_data_A[0].eig_cg_n_iterations = 100;
            }
          smoother_data_A[level].preconditioner = mg_matrices_A_block[level].get_matrix_diagonal_inverse();
        }
      mg_smoother_A.initialize(mg_matrices_A_block, smoother_data_A);
    }

    // Schur complement matrix GMG Smoother: Chebyshev, degree 4
    typedef PreconditionChebyshev<SchurComplementMatrixType,vector_t> MSmootherType;
    mg::SmootherRelaxation<MSmootherType, vector_t>
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
              smoother_data_Schur[0].degree = numbers::invalid_unsigned_int;
              smoother_data_Schur[0].eig_cg_n_iterations = 100; /*mg_matrices_M[0].m();*/
            }
          smoother_data_Schur[level].preconditioner = mg_matrices_Schur_complement[level].get_matrix_diagonal_inverse();
        }
      mg_smoother_Schur.initialize(mg_matrices_Schur_complement, smoother_data_Schur);
    }

    // Coarse Solver is just an application of the Chebyshev smoother setup
    // in such a way to be a solver
    //ABlock GMG
    MGCoarseGridApplySmoother<vector_t> mg_coarse_A;
    mg_coarse_A.initialize(mg_smoother_A);

    //Schur complement matrix GMG
    MGCoarseGridApplySmoother<vector_t> mg_coarse_Schur;
    mg_coarse_Schur.initialize(mg_smoother_Schur);

    // Interface matrices
    // Ablock GMG
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<ABlockMatrixType> > mg_interface_matrices_A;
    mg_interface_matrices_A.resize(0, sim.triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<sim.triangulation.n_global_levels(); ++level)
      mg_interface_matrices_A[level].initialize(mg_matrices_A_block[level]);
    mg::Matrix<vector_t > mg_interface_A(mg_interface_matrices_A);

    // Schur complement matrix GMG
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<SchurComplementMatrixType> > mg_interface_matrices_Schur;
    mg_interface_matrices_Schur.resize(0, sim.triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<sim.triangulation.n_global_levels(); ++level)
      mg_interface_matrices_Schur[level].initialize(mg_matrices_Schur_complement[level]);
    mg::Matrix<vector_t > mg_interface_Schur(mg_interface_matrices_Schur);

    // MG Matrix
    mg::Matrix<vector_t > mg_matrix_A(mg_matrices_A_block);
    mg::Matrix<vector_t > mg_matrix_Schur(mg_matrices_Schur_complement);

    // MG object
    // ABlock GMG
    Multigrid<vector_t > mg_A(mg_matrix_A,
                              mg_coarse_A,
                              mg_transfer_A_block,
                              mg_smoother_A,
                              mg_smoother_A);
    mg_A.set_edge_matrices(mg_interface_A, mg_interface_A);

    // Schur complement matrix GMG
    Multigrid<vector_t > mg_Schur(mg_matrix_Schur,
                                  mg_coarse_Schur,
                                  mg_transfer_Schur_complement,
                                  mg_smoother_Schur,
                                  mg_smoother_Schur);
    mg_Schur.set_edge_matrices(mg_interface_Schur, mg_interface_Schur);

    // GMG Preconditioner for ABlock and Schur complement
    typedef PreconditionMG<dim, vector_t, MGTransferMatrixFree<dim,double> > GMGPreconditioner;
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
                                  linearized_stokes_initial_guess,
                                  sim.current_linearization_point);
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

        solution_copy.collect_sizes();
        initial_copy.collect_sizes();
        rhs_copy.collect_sizes();

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

    solution_copy.collect_sizes();
    rhs_copy.collect_sizes();

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
                          sim.parameters.linear_solver_A_block_tolerance,
                          sim.parameters.linear_solver_S_block_tolerance);

    // create an expensive preconditioner that solves for the A block with CG
    const internal::BlockSchurGMGPreconditioner<StokesMatrixType, ABlockMatrixType, SchurComplementMatrixType, GMGPreconditioner, GMGPreconditioner>
    preconditioner_expensive (stokes_matrix, A_block_matrix, Schur_complement_block_matrix,
                              prec_A, prec_Schur,
                              /*do_solve_A*/true,
                              /*do_solve_Schur*/true,
                              sim.parameters.linear_solver_A_block_tolerance,
                              sim.parameters.linear_solver_S_block_tolerance);

    PrimitiveVectorMemory<dealii::LinearAlgebra::distributed::BlockVector<double> > mem;

    // step 1a: try if the simple and fast solver
    // succeeds in n_cheap_stokes_solver_steps steps or less.
    try
      {
        // if this cheaper solver is not desired, then simply short-cut
        // the attempt at solving with the cheaper preconditioner
        if (sim.parameters.n_cheap_stokes_solver_steps == 0)
          throw SolverControl::NoConvergence(0,0);

        // Unlike with the expensive preconditioner which uses CG solves on both the
        // velocity and pressure space, the cheap preonditioner only contains matrix-vector
        // products and GMG v-cycle where the smoothers, transfer operators, and coarse
        // solvers are all defined to be linear operators which do not change from iteration
        // to iterations. Therefore we can use GMRES, instead of FGMRES, as the Krylov subspace solver.
        SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double> >
        solver(solver_control_cheap, mem,
               SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double> >::
               AdditionalData(sim.parameters.stokes_gmres_restart_length+2,
                              true));

        solver.solve (stokes_matrix,
                      solution_copy,
                      rhs_copy,
                      preconditioner_cheap);

        final_linear_residual = solver_control_cheap.last_value();
      }
    // step 1b: take the stronger solver in case
    // the simple solver failed and attempt solving
    // it in n_expensive_stokes_solver_steps steps or less.
    catch (const SolverControl::NoConvergence &)
      {
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
            AssertThrow (sim.parameters.n_expensive_stokes_solver_steps>0,
                         ExcMessage ("The Stokes solver did not converge in the number of requested cheap iterations and "
                                     "you requested 0 for ``Maximum number of expensive Stokes solver steps''. Aborting."));

            solver.solve(stokes_matrix,
                         solution_copy,
                         rhs_copy,
                         preconditioner_expensive);

            final_linear_residual = solver_control_expensive.last_value();
          }
        // if the solver fails, report the error from processor 0 with some additional
        // information about its location, and throw a quiet exception on all other
        // processors
        catch (const std::exception &exc)
          {
            sim.signals.post_stokes_solver(sim,
                                           preconditioner_cheap.n_iterations_Schur_complement() + preconditioner_expensive.n_iterations_Schur_complement(),
                                           preconditioner_cheap.n_iterations_A_block() + preconditioner_expensive.n_iterations_A_block(),
                                           solver_control_cheap,
                                           solver_control_expensive);

            if (Utilities::MPI::this_mpi_process(sim.mpi_communicator) == 0)
              {
                // output solver history
                std::ofstream f((sim.parameters.output_directory+"solver_history.txt").c_str());

                // Only request the solver history if a history has actually been created
                if (sim.parameters.n_cheap_stokes_solver_steps > 0)
                  {
                    for (unsigned int i=0; i<solver_control_cheap.get_history_data().size(); ++i)
                      f << i << " " << solver_control_cheap.get_history_data()[i] << "\n";

                    f << "\n";
                  }


                for (unsigned int i=0; i<solver_control_expensive.get_history_data().size(); ++i)
                  f << i << " " << solver_control_expensive.get_history_data()[i] << "\n";

                f.close();

                AssertThrow (false,
                             ExcMessage (std::string("The iterative Stokes solver "
                                                     "did not converge. It reported the following error:\n\n")
                                         +
                                         exc.what()
                                         + "\n See " + sim.parameters.output_directory+"solver_history.txt"
                                         + " for convergence history."));
              }
            else
              {
                throw QuietException();
              }
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

    sim.current_constraints.distribute (distributed_stokes_solution);

    // now rescale the pressure back to real physical units
    distributed_stokes_solution.block(block_p) *= sim.pressure_scaling;

    // then copy back the solution from the temporary (non-ghosted) vector
    // into the ghosted one with all solution components
    sim.solution.block(block_vel) = distributed_stokes_solution.block(block_vel);
    sim.solution.block(block_p) = distributed_stokes_solution.block(block_p);

    // print the number of iterations to screen
    sim.pcout << (solver_control_cheap.last_step() != numbers::invalid_unsigned_int ?
                  solver_control_cheap.last_step():
                  0)
              << '+'
              << (solver_control_expensive.last_step() != numbers::invalid_unsigned_int ?
                  solver_control_expensive.last_step():
                  0)
              << " iterations.";
    sim.pcout << std::endl;

    // do some cleanup now that we have the solution
    sim.remove_nullspace(sim.solution, distributed_stokes_solution);
    if (sim.assemble_newton_stokes_system == false)
      sim.last_pressure_normalization_adjustment = sim.normalize_pressure(sim.solution);


    // convert melt pressures
    // TODO: We assert in the StokesMatrixFreeHandler constructor that we
    //       are not including melt transport.
    if (sim.parameters.include_melt_transport)
      sim.melt_handler->compute_melt_variables(sim.solution);

    return std::pair<double,double>(initial_nonlinear_residual,
                                    final_linear_residual);
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::setup_dofs()
  {
    // Velocity DoFHandler
    {
      dof_handler_v.clear();
      dof_handler_v.distribute_dofs(fe_v);

      DoFRenumbering::hierarchical(dof_handler_v);

      constraints_v.clear();
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs (dof_handler_v,
                                               locally_relevant_dofs);
      constraints_v.reinit(locally_relevant_dofs);
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

      constraints_p.clear();
      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs (dof_handler_p,
                                               locally_relevant_dofs);
      constraints_p.reinit(locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints (dof_handler_p, constraints_p);
      constraints_p.close();
    }

    // Coefficient transfer objects
    {
      dof_handler_projection.clear();
      dof_handler_projection.distribute_dofs(fe_projection);

      DoFRenumbering::hierarchical(dof_handler_projection);

      active_coef_dof_vec.reinit(dof_handler_projection.locally_owned_dofs(), sim.triangulation.get_communicator());
    }

    // Multigrid DoF setup
    {
      //Ablock GMG
      dof_handler_v.distribute_mg_dofs();

      mg_constrained_dofs_A_block.clear();
      mg_constrained_dofs_A_block.initialize(dof_handler_v);

      std::set<types::boundary_id> dirichlet_boundary = sim.boundary_velocity_manager.get_zero_boundary_velocity_indicators();
      for (auto it: sim.boundary_velocity_manager.get_active_boundary_velocity_names())
        {
          types::boundary_id bdryid = it.first;
          std::string component=it.second.first;
          Assert(component=="", ExcNotImplemented());
          dirichlet_boundary.insert(bdryid);
        }
      mg_constrained_dofs_A_block.make_zero_boundary_constraints(dof_handler_v, dirichlet_boundary);

      {
        std::set<types::boundary_id> no_flux_boundary = sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
        if (!no_flux_boundary.empty() && !sim.geometry_model->has_curved_elements())
          for (auto bid : no_flux_boundary)
            {
              internal::TangentialBoundaryFunctions::compute_no_normal_flux_constraints_box(dof_handler_v,
                                                                                            bid,
                                                                                            0,
                                                                                            mg_constrained_dofs_A_block);
            }
      }

      //Schur complement matrix GMG
      dof_handler_p.distribute_mg_dofs();

      mg_constrained_dofs_Schur_complement.clear();
      mg_constrained_dofs_Schur_complement.initialize(dof_handler_p);

      dof_handler_projection.distribute_mg_dofs();
    }

    // Setup the matrix-free operators
    // Stokes matrix
    {
      typename MatrixFree<dim,double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim,double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_gradients |
                                              update_JxW_values | update_quadrature_points);

      std::vector<const DoFHandler<dim>*> stokes_dofs;
      stokes_dofs.push_back(&dof_handler_v);
      stokes_dofs.push_back(&dof_handler_p);
      std::vector<const ConstraintMatrix *> stokes_constraints;
      stokes_constraints.push_back(&constraints_v);
      stokes_constraints.push_back(&constraints_p);

      std::shared_ptr<MatrixFree<dim,double> >
      stokes_mf_storage(new MatrixFree<dim,double>());
      stokes_mf_storage->reinit(*sim.mapping,stokes_dofs, stokes_constraints,
                                QGauss<1>(sim.parameters.stokes_velocity_degree+1), additional_data);
      stokes_matrix.clear();
      stokes_matrix.initialize(stokes_mf_storage);

    }

    // ABlock matrix
    {
      typename MatrixFree<dim,double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim,double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_gradients |
                                              update_JxW_values | update_quadrature_points);
      std::shared_ptr<MatrixFree<dim,double> >
      ablock_mf_storage(new MatrixFree<dim,double>());
      ablock_mf_storage->reinit(*sim.mapping,dof_handler_v, constraints_v,
                                QGauss<1>(sim.parameters.stokes_velocity_degree+1), additional_data);

      A_block_matrix.clear();
      A_block_matrix.initialize(ablock_mf_storage);
    }

    // Schur complement block matrix
    {
      typename MatrixFree<dim,double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim,double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_values | update_JxW_values |
                                              update_quadrature_points);
      std::shared_ptr<MatrixFree<dim,double> >
      Schur_mf_storage(new MatrixFree<dim,double>());
      Schur_mf_storage->reinit(*sim.mapping,dof_handler_p, constraints_p,
                               QGauss<1>(sim.parameters.stokes_velocity_degree+1), additional_data);

      Schur_complement_block_matrix.clear();
      Schur_complement_block_matrix.initialize(Schur_mf_storage);
    }

    // GMG matrices
    {
      const unsigned int n_levels = sim.triangulation.n_global_levels();

      // ABlock GMG
      mg_matrices_A_block.clear_elements();
      mg_matrices_A_block.resize(0, n_levels-1);

      for (unsigned int level=0; level<n_levels; ++level)
        {
          IndexSet relevant_dofs;
          DoFTools::extract_locally_relevant_level_dofs(dof_handler_v, level, relevant_dofs);
          ConstraintMatrix level_constraints;
          level_constraints.reinit(relevant_dofs);
          level_constraints.add_lines(mg_constrained_dofs_A_block.get_boundary_indices(level));
          level_constraints.close();

          std::set<types::boundary_id> no_flux_boundary
            = sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
          if (!no_flux_boundary.empty() && sim.geometry_model->has_curved_elements())
            {
#if DEAL_II_VERSION_GTE(9,2,0)
              ConstraintMatrix user_level_constraints;
              user_level_constraints.reinit(relevant_dofs);

              internal::TangentialBoundaryFunctions::compute_no_normal_flux_constraints_shell(dof_handler_v,
                                                                                              mg_constrained_dofs_A_block,
                                                                                              *sim.mapping,
                                                                                              level,
                                                                                              0,
                                                                                              no_flux_boundary,
                                                                                              user_level_constraints);
              user_level_constraints.close();
              mg_constrained_dofs_A_block.add_user_constraints(level,user_level_constraints);

              // let Dirichlet values win over no normal flux:
              level_constraints.merge(user_level_constraints, ConstraintMatrix::left_object_wins);
              level_constraints.close();
#else
              AssertThrow(false, ExcMessage("No normal flux for spherical domains requires "
                                            "a deal.II version newer than 9.1"));
#endif
            }

          {
            typename MatrixFree<dim,double>::AdditionalData additional_data;
            additional_data.tasks_parallel_scheme =
              MatrixFree<dim,double>::AdditionalData::none;
            additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                                    update_quadrature_points);
#if DEAL_II_VERSION_GTE(9,2,0)
            additional_data.mg_level = level;
#else
            additional_data.level_mg_handler = level;
#endif
            std::shared_ptr<MatrixFree<dim,double> >
            mg_mf_storage_level(new MatrixFree<dim,double>());
            mg_mf_storage_level->reinit(*sim.mapping, dof_handler_v, level_constraints,
                                        QGauss<1>(sim.parameters.stokes_velocity_degree+1),
                                        additional_data);

            mg_matrices_A_block[level].clear();
            mg_matrices_A_block[level].initialize(mg_mf_storage_level, mg_constrained_dofs_A_block, level);

          }
        }

      //Schur complement matrix GMG
      mg_matrices_Schur_complement.clear_elements();
      mg_matrices_Schur_complement.resize(0, n_levels-1);

      for (unsigned int level=0; level<n_levels; ++level)
        {
          IndexSet relevant_dofs;
          DoFTools::extract_locally_relevant_level_dofs(dof_handler_p, level, relevant_dofs);
          ConstraintMatrix level_constraints;
          level_constraints.reinit(relevant_dofs);
          level_constraints.close();

          {
            typename MatrixFree<dim,double>::AdditionalData additional_data;
            additional_data.tasks_parallel_scheme =
              MatrixFree<dim,double>::AdditionalData::none;
            additional_data.mapping_update_flags = (update_values | update_JxW_values |
                                                    update_quadrature_points);
#if DEAL_II_VERSION_GTE(9,2,0)
            additional_data.mg_level = level;
#else
            additional_data.level_mg_handler = level;
#endif
            std::shared_ptr<MatrixFree<dim,double> >
            mg_mf_storage_level(new MatrixFree<dim,double>());
            mg_mf_storage_level->reinit(*sim.mapping, dof_handler_p, level_constraints,
                                        QGauss<1>(sim.parameters.stokes_velocity_degree+1),
                                        additional_data);

            mg_matrices_Schur_complement[level].clear();
            mg_matrices_Schur_complement[level].initialize(mg_mf_storage_level, mg_constrained_dofs_Schur_complement, level);
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

    // GMG diagonals
    for (unsigned int level=0; level < sim.triangulation.n_global_levels(); ++level)
      {
        mg_matrices_Schur_complement[level].compute_diagonal();

        // If we have a tangential boundary we must compute the A block
        // diagonal outside of the matrix-free object
        if (!(sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators().empty())
            &&
            sim.geometry_model->has_curved_elements())
          {
            IndexSet locally_relevant_dofs;
            DoFTools::extract_locally_relevant_level_dofs (dof_handler_v, level, locally_relevant_dofs);

            DiagonalMatrix<dealii::LinearAlgebra::distributed::Vector<double> > diagonal_matrix;
            dealii::LinearAlgebra::distributed::Vector<double> &diagonal_vector =
              diagonal_matrix.get_vector();

            diagonal_vector.reinit(dof_handler_v.locally_owned_mg_dofs(level),
                                   locally_relevant_dofs,
                                   sim.mpi_communicator);

            QGauss<dim>  quadrature_formula(sim.parameters.stokes_velocity_degree+1);
            FEValues<dim> fe_values (fe_v, quadrature_formula,
                                     update_values   | update_gradients |
                                     update_quadrature_points | update_JxW_values);

            const unsigned int   dofs_per_cell   = fe_v.dofs_per_cell;
            const unsigned int   n_q_points      = quadrature_formula.size();

            FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

            std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
            const FEValuesExtractors::Vector velocities (0);

            std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);

            ConstraintMatrix boundary_constraints;
            boundary_constraints.reinit(locally_relevant_dofs);
            boundary_constraints.add_lines (mg_constrained_dofs_A_block.get_refinement_edge_indices(level));
            boundary_constraints.add_lines (mg_constrained_dofs_A_block.get_boundary_indices(level));
#if DEAL_II_VERSION_GTE(9,2,0)
            // let Dirichlet values win over no normal flux:
            boundary_constraints.merge(mg_constrained_dofs_A_block.get_user_constraint_matrix(level),
                                       ConstraintMatrix::left_object_wins);
#endif
            boundary_constraints.close();

            typename DoFHandler<dim>::level_cell_iterator cell = dof_handler_v.begin(level),
                                                          endc = dof_handler_v.end(level);
            for (; cell!=endc; ++cell)
              if (cell->level_subdomain_id()==sim.triangulation.locally_owned_subdomain())
                {
                  cell_matrix = 0;
                  fe_values.reinit (cell);

                  typename DoFHandler<dim>::level_cell_iterator DG_cell(&(sim.triangulation),
                                                                        level,
                                                                        cell->index(),
                                                                        &dof_handler_projection);
                  std::vector<types::global_dof_index> dg_dof_indices(dof_handler_projection.get_fe(0).dofs_per_cell);
                  DG_cell->get_active_or_mg_dof_indices(dg_dof_indices);
                  double viscosity = level_coef_dof_vec[level](dg_dof_indices[0]);

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      for (unsigned int k=0; k<dofs_per_cell; ++k)
                        symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);

                      const double JxW = fe_values.JxW(q);
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          cell_matrix(i,j) += 2. * viscosity * (symgrad_phi_u[i]*symgrad_phi_u[j])
                                              * JxW;
                    }

                  cell->get_mg_dof_indices (local_dof_indices);

                  boundary_constraints.distribute_local_to_global (cell_matrix,
                                                                   local_dof_indices,
                                                                   diagonal_matrix);
                }

            mg_matrices_A_block[level].set_diagonal(diagonal_matrix.get_vector());
          }
        else
          {
            mg_matrices_A_block[level].compute_diagonal();
          }
      }
  }




// explicit instantiation of the functions we implement in this file
#define INSTANTIATE(dim) \
  template class StokesMatrixFreeHandler<dim>; \
  template class StokesMatrixFreeHandlerImplementation<dim,2>; \
  template class StokesMatrixFreeHandlerImplementation<dim,3>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
