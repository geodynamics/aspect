/*
  Copyright (C) 2018 - 2020 by the authors of the ASPECT code.

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

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/read_write_vector.templates.h>

#include <deal.II/lac/solver_idr.h>



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
                     AffineConstraints<double> &constraints,
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
                                                    AffineConstraints<double> &constraints)
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
    viscosity = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::BlockVector<number> >::clear();
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::StokesOperator<dim,degree_v,number>::
  fill_cell_data (const Table<2, VectorizedArray<number>> &viscosity_table,
                  const double pressure_scaling,
                  const bool is_compressible)
  {
    viscosity = &viscosity_table;
    this->pressure_scaling = pressure_scaling;
    this->is_compressible = is_compressible;
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

    const bool use_viscosity_at_quadrature_points
      = (viscosity->size(1) == velocity.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        VectorizedArray<number> viscosity_x_2 = 2.0*(*viscosity)(cell, 0);

        velocity.reinit (cell);
        velocity.read_dof_values (src.block(0));
        velocity.evaluate (false,true,false);
        pressure.reinit (cell);
        pressure.read_dof_values (src.block(1));
        pressure.evaluate (true,false,false);

        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              viscosity_x_2 = 2.0*(*viscosity)(cell, q);

            SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
                                                          velocity.get_symmetric_gradient (q);
            VectorizedArray<number> pres = pressure.get_value(q);
            VectorizedArray<number> div = trace(sym_grad_u);
            pressure.submit_value(-pressure_scaling*div, q);

            sym_grad_u *= viscosity_x_2;

            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= pressure_scaling*pres;

            if (is_compressible)
              for (unsigned int d=0; d<dim; ++d)
                sym_grad_u[d][d] -= viscosity_x_2/3.0*div;

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
    viscosity = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::clear();
  }

  template <int dim, int degree_p, typename number>
  void
  MatrixFreeStokesOperators::MassMatrixOperator<dim,degree_p,number>::
  fill_cell_data (const Table<2, VectorizedArray<number>> &viscosity_table,
                  const double pressure_scaling)
  {
    viscosity = &viscosity_table;
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

    const bool use_viscosity_at_quadrature_points
      = (viscosity->size(1) == pressure.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        VectorizedArray<number> one_over_viscosity = (*viscosity)(cell, 0);

#if DEAL_II_VERSION_GTE(9,3,0)
        const unsigned int n_components_filled = this->get_matrix_free()->n_active_entries_per_cell_batch(cell);
#else
        const unsigned int n_components_filled = this->get_matrix_free()->n_components_filled(cell);
#endif

        // The /= operator for VectorizedArray results in a floating point operation
        // (divide by 0) since the (*viscosity)(cell) array is not completely filled.
        // Therefore, we need to divide each entry manually.
        for (unsigned int c=0; c<n_components_filled; ++c)
          one_over_viscosity[c] = pressure_scaling*pressure_scaling/one_over_viscosity[c];

        pressure.reinit (cell);
        pressure.read_dof_values(src);
        pressure.evaluate (true, false);
        for (unsigned int q=0; q<pressure.n_q_points; ++q)
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              {
                one_over_viscosity = (*viscosity)(cell, q);

#if DEAL_II_VERSION_GTE(9,3,0)
                const unsigned int n_components_filled = this->get_matrix_free()->n_active_entries_per_cell_batch(cell);
#else
                const unsigned int n_components_filled = this->get_matrix_free()->n_components_filled(cell);
#endif

                for (unsigned int c=0; c<n_components_filled; ++c)
                  one_over_viscosity[c] = pressure_scaling*pressure_scaling/one_over_viscosity[c];
              }

            pressure.submit_value(one_over_viscosity*
                                  pressure.get_value(q),q);
          }
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

    const bool use_viscosity_at_quadrature_points
      = (viscosity->size(1) == pressure.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        VectorizedArray<number> one_over_viscosity = (*viscosity)(cell, 0);

#if DEAL_II_VERSION_GTE(9,3,0)
        const unsigned int n_components_filled = this->get_matrix_free()->n_active_entries_per_cell_batch(cell);
#else
        const unsigned int n_components_filled = this->get_matrix_free()->n_components_filled(cell);
#endif

        // The /= operator for VectorizedArray results in a floating point operation
        // (divide by 0) since the (*viscosity)(cell) array is not completely filled.
        // Therefore, we need to divide each entry manually.
        for (unsigned int c=0; c<n_components_filled; ++c)
          one_over_viscosity[c] = pressure_scaling*pressure_scaling/one_over_viscosity[c];

        pressure.reinit (cell);
        AlignedVector<VectorizedArray<number> > diagonal(pressure.dofs_per_cell);
        for (unsigned int i=0; i<pressure.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<pressure.dofs_per_cell; ++j)
              pressure.begin_dof_values()[j] = VectorizedArray<number>();
            pressure.begin_dof_values()[i] = make_vectorized_array<number> (1.);

            pressure.evaluate (true,false,false);
            for (unsigned int q=0; q<pressure.n_q_points; ++q)
              {
                // Only update the viscosity if a Q1 projection is used.
                if (use_viscosity_at_quadrature_points)
                  {
                    one_over_viscosity = (*viscosity)(cell, q);

#if DEAL_II_VERSION_GTE(9,3,0)
                    const unsigned int n_components_filled = this->get_matrix_free()->n_active_entries_per_cell_batch(cell);
#else
                    const unsigned int n_components_filled = this->get_matrix_free()->n_components_filled(cell);
#endif

                    for (unsigned int c=0; c<n_components_filled; ++c)
                      one_over_viscosity[c] = pressure_scaling*pressure_scaling/one_over_viscosity[c];
                  }

                pressure.submit_value(one_over_viscosity*
                                      pressure.get_value(q),q);
              }
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
    viscosity = nullptr;
    MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::clear();
  }

  template <int dim, int degree_v, typename number>
  void
  MatrixFreeStokesOperators::ABlockOperator<dim,degree_v,number>::
  fill_cell_data (const Table<2, VectorizedArray<number>> &viscosity_table,
                  const bool is_compressible)
  {
    viscosity = &viscosity_table;
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

    const bool use_viscosity_at_quadrature_points
      = (viscosity->size(1) == velocity.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        VectorizedArray<number> viscosity_x_2 = 2.0*(*viscosity)(cell, 0);

        velocity.reinit (cell);
        velocity.read_dof_values(src);
        velocity.evaluate (false, true, false);
        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              viscosity_x_2 = 2.0*(*viscosity)(cell, q);

            SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
                                                          velocity.get_symmetric_gradient (q);
            sym_grad_u *= viscosity_x_2;

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

    const bool use_viscosity_at_quadrature_points
      = (viscosity->size(1) == velocity.n_q_points);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        VectorizedArray<number> viscosity_x_2 = 2.0*(*viscosity)(cell, 0);

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
                // Only update the viscosity if a Q1 projection is used.
                if (use_viscosity_at_quadrature_points)
                  viscosity_x_2 = 2.0*(*viscosity)(cell, q);

                SymmetricTensor<2,dim,VectorizedArray<number>> sym_grad_u =
                                                              velocity.get_symmetric_gradient (q);

                sym_grad_u *= viscosity_x_2;

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
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::declare_parameters(ParameterHandler &/*prm*/)
  {
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim,velocity_degree>::parse_parameters(ParameterHandler &/*prm*/)
  {
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
                  | avg::pick_largest | avg::project_to_Q1 | avg::log_average
                  | avg::harmonic_average_only_viscosity | avg::project_to_Q1_only_viscosity)) != 0,
                ExcMessage("The matrix-free Stokes solver currently only works if material model averaging "
                           "is enabled. If no averaging is desired, consider using ``project to Q1 only "
                           "viscosity''."));

    // Currently cannot solve compressible flow with implicit reference density
    if (sim.material_model->is_compressible() == true)
      AssertThrow(sim.parameters.formulation_mass_conservation !=
                  Parameters<dim>::Formulation::MassConservation::implicit_reference_density_profile,
                  ExcNotImplemented());

    {
      const unsigned int n_vect_doubles =
        VectorizedArray<double>::size();
      const unsigned int n_vect_bits = 8 * sizeof(double) * n_vect_doubles;

      sim.pcout << "Vectorization over " << n_vect_doubles
                << " doubles = " << n_vect_bits << " bits ("
                << dealii::Utilities::System::get_current_vectorization_level()
                << "), VECTORIZATION_LEVEL=" << DEAL_II_COMPILER_VECTORIZATION_LEVEL
                << std::endl;
    }
  }


  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::evaluate_material_model ()
  {
    dealii::LinearAlgebra::distributed::Vector<double> active_viscosity_vector(dof_handler_projection.locally_owned_dofs(),
                                                                               sim.triangulation.get_communicator());

    const QGauss<dim> quadrature_formula (sim.parameters.stokes_velocity_degree+1);

    double min_el = std::numeric_limits<double>::max();
    double max_el = -std::numeric_limits<double>::max();

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
          out);

        for (unsigned int i=0; i<values.size(); ++i)
          {
            // Find the max/min of the evaluated viscosities.
            min_el = std::min(min_el, out.viscosities[i]);
            max_el = std::max(max_el, out.viscosities[i]);

            values[i] = out.viscosities[i];
          }
        return;
      },
      active_viscosity_vector);

      active_viscosity_vector.compress(VectorOperation::insert);
    }

    FEValues<dim> fe_values_projection (*(sim.mapping),
                                        fe_projection,
                                        quadrature_formula,
                                        update_values);

    // Create active mesh viscosity table.
    {
#if DEAL_II_VERSION_GTE(9,3,0)
      const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_cell_batches();
#else
      const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_macro_cells();
#endif

      const unsigned int n_q_points = quadrature_formula.size();

      std::vector<double> values_on_quad;

      // One value per cell is required for DGQ0 projection and n_q_points
      // values per cell for DGQ1.
      if (dof_handler_projection.get_fe().degree == 0)
        active_viscosity_table.reinit(TableIndices<2>(n_cells, 1));
      else if (dof_handler_projection.get_fe().degree == 1)
        {
          values_on_quad.resize(n_q_points);
          active_viscosity_table.reinit(TableIndices<2>(n_cells, n_q_points));
        }
      else
        Assert(false, ExcInternalError());

      std::vector<types::global_dof_index> local_dof_indices(fe_projection.dofs_per_cell);
      for (unsigned int cell=0; cell<n_cells; ++cell)
        {
#if DEAL_II_VERSION_GTE(9,3,0)
          const unsigned int n_components_filled = stokes_matrix.get_matrix_free()->n_active_entries_per_cell_batch(cell);
#else
          const unsigned int n_components_filled = stokes_matrix.get_matrix_free()->n_components_filled(cell);
#endif

          for (unsigned int i=0; i<n_components_filled; ++i)
            {
              typename DoFHandler<dim>::active_cell_iterator FEQ_cell =
                stokes_matrix.get_matrix_free()->get_cell_iterator(cell,i);
              typename DoFHandler<dim>::active_cell_iterator DG_cell(&(sim.triangulation),
                                                                     FEQ_cell->level(),
                                                                     FEQ_cell->index(),
                                                                     &dof_handler_projection);
              DG_cell->get_active_or_mg_dof_indices(local_dof_indices);

              // For DGQ0, we simply use the viscosity at the single
              // support point of the element. For DGQ1, we must project
              // back to quadrature point values.
              if (dof_handler_projection.get_fe().degree == 0)
                active_viscosity_table(cell, 0)[i] = active_viscosity_vector(local_dof_indices[0]);
              else
                {
                  fe_values_projection.reinit(DG_cell);
                  fe_values_projection.get_function_values(active_viscosity_vector,
                                                           local_dof_indices,
                                                           values_on_quad);

                  // Do not allow viscosity to be greater than or less than the limits
                  // of the evaluated viscosity on the active level.
                  for (unsigned int q=0; q<n_q_points; ++q)
                    active_viscosity_table(cell, q)[i]
                      = std::min(std::max(values_on_quad[q], min_el), max_el);
                }
            }
        }
    }

    const bool is_compressible = sim.material_model->is_compressible();

    // Store viscosity tables and other data into the active level matrix-free objects.
    stokes_matrix.fill_cell_data(active_viscosity_table,
                                 sim.pressure_scaling,
                                 is_compressible);

    if (sim.parameters.n_expensive_stokes_solver_steps > 0)
      {
        A_block_matrix.fill_cell_data(active_viscosity_table,
                                      is_compressible);
        Schur_complement_block_matrix.fill_cell_data(active_viscosity_table,
                                                     sim.pressure_scaling);
      }

    const unsigned int n_levels = sim.triangulation.n_global_levels();
    level_viscosity_vector = 0.;
    level_viscosity_vector.resize(0,n_levels-1);

    // Project the active level viscosity vector to multilevel vector representations
    // using MG transfer objects. This transfer is based on the same linear operator used to
    // transfer data inside a v-cycle.
    MGTransferMatrixFree<dim,GMGNumberType> transfer;
    transfer.build(dof_handler_projection);

    // Explicitly pick the version with template argument double to convert
    // double-valued active_viscosity_vector to GMGNumberType-valued
    // level_viscosity_vector:
    transfer.template interpolate_to_mg<double>(dof_handler_projection,
                                                level_viscosity_vector,
                                                active_viscosity_vector);

    level_viscosity_tables.resize(0,n_levels-1);
    for (unsigned int level=0; level<n_levels; ++level)
      {
        // Create viscosity tables on each level.
#if DEAL_II_VERSION_GTE(9,3,0)
        const unsigned int n_cells = mg_matrices_A_block[level].get_matrix_free()->n_cell_batches();
#else
        const unsigned int n_cells = mg_matrices_A_block[level].get_matrix_free()->n_macro_cells();
#endif

        const unsigned int n_q_points = quadrature_formula.size();

        std::vector<GMGNumberType> values_on_quad;

        // One value per cell is required for DGQ0 projection and n_q_points
        // values per cell for DGQ1.
        if (dof_handler_projection.get_fe().degree == 0)
          level_viscosity_tables[level].reinit(TableIndices<2>(n_cells, 1));
        else
          {
            values_on_quad.resize(n_q_points);
            level_viscosity_tables[level].reinit(TableIndices<2>(n_cells, n_q_points));
          }

        std::vector<types::global_dof_index> local_dof_indices(fe_projection.dofs_per_cell);
        for (unsigned int cell=0; cell<n_cells; ++cell)
          {
#if DEAL_II_VERSION_GTE(9,3,0)
            const unsigned int n_components_filled = mg_matrices_A_block[level].get_matrix_free()->n_active_entries_per_cell_batch(cell);
#else
            const unsigned int n_components_filled = mg_matrices_A_block[level].get_matrix_free()->n_components_filled(cell);
#endif

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
                  level_viscosity_tables[level](cell, 0)[i] = level_viscosity_vector[level](local_dof_indices[0]);
                else
                  {
                    fe_values_projection.reinit(DG_cell);
                    fe_values_projection.get_function_values(level_viscosity_vector[level],
                                                             local_dof_indices,
                                                             values_on_quad);

                    // Do not allow viscosity to be greater than or less than the limits
                    // of the evaluated viscosity on the active level.
                    for (unsigned int q=0; q<n_q_points; ++q)
                      level_viscosity_tables[level](cell,q)[i]
                        = std::min(std::max(values_on_quad[q], static_cast<GMGNumberType>(min_el)),
                                   static_cast<GMGNumberType>(max_el));
                  }
              }
          }

        // Store viscosity tables and other data into the multigrid level matrix-free objects.
        mg_matrices_A_block[level].fill_cell_data (level_viscosity_tables[level],
                                                   is_compressible);
        mg_matrices_Schur_complement[level].fill_cell_data (level_viscosity_tables[level],
                                                            sim.pressure_scaling);
      }
  }



  template <int dim, int velocity_degree>
  void StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::correct_stokes_rhs()
  {
    const bool is_compressible = sim.material_model->is_compressible();

    dealii::LinearAlgebra::distributed::BlockVector<double> rhs_correction(2);
    dealii::LinearAlgebra::distributed::BlockVector<double> u0(2);

    stokes_matrix.initialize_dof_vector(rhs_correction);
    stokes_matrix.initialize_dof_vector(u0);

    // The vector u0 is a zero vector, but with correct boundary values.
    u0 = 0;
    rhs_correction = 0;
    sim.current_constraints.distribute(u0);
    u0.update_ghost_values();

    FEEvaluation<dim,velocity_degree,velocity_degree+1,dim,double>
    velocity (*stokes_matrix.get_matrix_free(), 0);
    FEEvaluation<dim,velocity_degree-1,velocity_degree+1,1,double>
    pressure (*stokes_matrix.get_matrix_free(), 1);

    const bool use_viscosity_at_quadrature_points
      = (active_viscosity_table.size(1) == velocity.n_q_points);

#if DEAL_II_VERSION_GTE(9,3,0)
    const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_cell_batches();
#else
    const unsigned int n_cells = stokes_matrix.get_matrix_free()->n_macro_cells();
#endif

    // Much like the matrix-free apply_add() functions compute a matrix-vector
    // product by looping over cells and applying local matrix operations,
    // here we apply the negative of the stokes_matrix operator to u0.
    for (unsigned int cell=0; cell<n_cells; ++cell)
      {
        VectorizedArray<double> viscosity_x_2 = 2.0*active_viscosity_table(cell, 0);

        // We must use read_dof_values_plain() as to not overwrite boundary information
        // with the zero boundary used by the stokes_matrix operator.
        velocity.reinit (cell);
        velocity.read_dof_values_plain (u0.block(0));
        velocity.evaluate (false,true,false);
        pressure.reinit (cell);
        pressure.read_dof_values_plain (u0.block(1));
        pressure.evaluate (true,false,false);

        for (unsigned int q=0; q<velocity.n_q_points; ++q)
          {
            // Only update the viscosity if a Q1 projection is used.
            if (use_viscosity_at_quadrature_points)
              viscosity_x_2 = 2.0*active_viscosity_table(cell, q);

            SymmetricTensor<2,dim,VectorizedArray<double>> sym_grad_u =
                                                          velocity.get_symmetric_gradient (q);
            VectorizedArray<double> pres = pressure.get_value(q);
            VectorizedArray<double> div = trace(sym_grad_u);
            pressure.submit_value   (sim.pressure_scaling*div, q);

            sym_grad_u *= viscosity_x_2;

            for (unsigned int d=0; d<dim; ++d)
              sym_grad_u[d][d] -= sim.pressure_scaling*pres;

            if (is_compressible)
              for (unsigned int d=0; d<dim; ++d)
                sym_grad_u[d][d] -= viscosity_x_2/3.0*div;

            velocity.submit_symmetric_gradient(-1.0*sym_grad_u, q);
          }

        velocity.integrate (false,true);
        velocity.distribute_local_to_global (rhs_correction.block(0));
        pressure.integrate (true,false);
        pressure.distribute_local_to_global (rhs_correction.block(1));
      }
    rhs_correction.compress(VectorOperation::add);

    // Copy to the correct vector type and add the correction to the system rhs.
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
      }


    // Coarse Solver is just an application of the Chebyshev smoother setup
    // in such a way to be a solver
    //ABlock GMG
    MGCoarseGridApplySmoother<VectorType> mg_coarse_A;
    mg_coarse_A.initialize(mg_smoother_A);

    //Schur complement matrix GMG
    MGCoarseGridApplySmoother<VectorType> mg_coarse_Schur;
    mg_coarse_Schur.initialize(mg_smoother_Schur);

    // Interface matrices
    // Ablock GMG
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<GMGABlockMatrixType> > mg_interface_matrices_A;
    mg_interface_matrices_A.resize(0, sim.triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<sim.triangulation.n_global_levels(); ++level)
      mg_interface_matrices_A[level].initialize(mg_matrices_A_block[level]);
    mg::Matrix<VectorType > mg_interface_A(mg_interface_matrices_A);

    // Schur complement matrix GMG
    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<GMGSchurComplementMatrixType> > mg_interface_matrices_Schur;
    mg_interface_matrices_Schur.resize(0, sim.triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<sim.triangulation.n_global_levels(); ++level)
      mg_interface_matrices_Schur[level].initialize(mg_matrices_Schur_complement[level]);
    mg::Matrix<VectorType > mg_interface_Schur(mg_interface_matrices_Schur);

    // MG Matrix
    mg::Matrix<VectorType > mg_matrix_A(mg_matrices_A_block);
    mg::Matrix<VectorType > mg_matrix_Schur(mg_matrices_Schur_complement);

    // MG object
    // ABlock GMG
    Multigrid<VectorType > mg_A(mg_matrix_A,
                                mg_coarse_A,
                                mg_transfer_A_block,
                                mg_smoother_A,
                                mg_smoother_A);
    mg_A.set_edge_matrices(mg_interface_A, mg_interface_A);

    // Schur complement matrix GMG
    Multigrid<VectorType > mg_Schur(mg_matrix_Schur,
                                    mg_coarse_Schur,
                                    mg_transfer_Schur_complement,
                                    mg_smoother_Schur,
                                    mg_smoother_Schur);
    mg_Schur.set_edge_matrices(mg_interface_Schur, mg_interface_Schur);

    // GMG Preconditioner for ABlock and Schur complement
    using GMGPreconditioner = PreconditionMG<dim, VectorType, MGTransferMatrixFree<dim,GMGNumberType> >;
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
        // velocity and pressure space, the cheap preconditioner only contains matrix-vector
        // products and GMG v-cycle where the smoothers, transfer operators, and coarse
        // solvers are all defined to be linear operators which do not change from iteration
        // to iteration. Therefore we can use non-flexible Krylov methods like GMRES or IDR(s),
        // instead of requiring FGMRES, greatly lowing the memory requirement of the solver.
        if (sim.parameters.stokes_krylov_type == Parameters<dim>::StokesKrylovType::gmres)
          {
            SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double> >
            solver(solver_control_cheap, mem,
                   SolverGMRES<dealii::LinearAlgebra::distributed::BlockVector<double> >::
                   AdditionalData(sim.parameters.stokes_gmres_restart_length+2,
                                  true));

            solver.solve (stokes_matrix,
                          solution_copy,
                          rhs_copy,
                          preconditioner_cheap);
          }
        else if (sim.parameters.stokes_krylov_type == Parameters<dim>::StokesKrylovType::idr_s)
          {
            SolverIDR<dealii::LinearAlgebra::distributed::BlockVector<double> >
            solver(solver_control_cheap, mem,
                   SolverIDR<dealii::LinearAlgebra::distributed::BlockVector<double> >::
                   AdditionalData(sim.parameters.idr_s_parameter));

            solver.solve (stokes_matrix,
                          solution_copy,
                          rhs_copy,
                          preconditioner_cheap);
          }
        else
          Assert(false,ExcNotImplemented());

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
      sim.melt_handler->compute_melt_variables(sim.system_matrix,sim.solution,sim.system_rhs);


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
    }

    // Multigrid DoF setup
    {
      //Ablock GMG
      dof_handler_v.distribute_mg_dofs();

      mg_constrained_dofs_A_block.clear();
      mg_constrained_dofs_A_block.initialize(dof_handler_v);

      std::set<types::boundary_id> dirichlet_boundary = sim.boundary_velocity_manager.get_zero_boundary_velocity_indicators();
      for (const auto &it: sim.boundary_velocity_manager.get_active_boundary_velocity_names())
        {
          const types::boundary_id bdryid = it.first;
          const std::string component=it.second.first;
          Assert(component=="", ExcNotImplemented());
          dirichlet_boundary.insert(bdryid);
        }
      mg_constrained_dofs_A_block.make_zero_boundary_constraints(dof_handler_v, dirichlet_boundary);

      {
        const std::set<types::boundary_id> no_flux_boundary = sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
        if (!no_flux_boundary.empty() && !sim.geometry_model->has_curved_elements())
          for (const auto bid : no_flux_boundary)
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
      std::vector<const AffineConstraints<double> *> stokes_constraints;
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
          AffineConstraints<double> level_constraints;
          level_constraints.reinit(relevant_dofs);
          level_constraints.add_lines(mg_constrained_dofs_A_block.get_boundary_indices(level));
          level_constraints.close();

          std::set<types::boundary_id> no_flux_boundary
            = sim.boundary_velocity_manager.get_tangential_boundary_velocity_indicators();
          if (!no_flux_boundary.empty() && sim.geometry_model->has_curved_elements())
            {
              AffineConstraints<double> user_level_constraints;
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
              level_constraints.merge(user_level_constraints, AffineConstraints<double>::left_object_wins);
              level_constraints.close();
            }

          {
            typename MatrixFree<dim,GMGNumberType>::AdditionalData additional_data;
            additional_data.tasks_parallel_scheme =
              MatrixFree<dim,GMGNumberType>::AdditionalData::none;
            additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                                    update_quadrature_points);
            additional_data.mg_level = level;
            std::shared_ptr<MatrixFree<dim,GMGNumberType> >
            mg_mf_storage_level(new MatrixFree<dim,GMGNumberType>());
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
          AffineConstraints<double> level_constraints;
          level_constraints.reinit(relevant_dofs);
          level_constraints.close();

          {
            typename MatrixFree<dim,GMGNumberType>::AdditionalData additional_data;
            additional_data.tasks_parallel_scheme =
              MatrixFree<dim,GMGNumberType>::AdditionalData::none;
            additional_data.mapping_update_flags = (update_values | update_JxW_values |
                                                    update_quadrature_points);
            additional_data.mg_level = level;
            std::shared_ptr<MatrixFree<dim,GMGNumberType> >
            mg_mf_storage_level(new MatrixFree<dim,GMGNumberType>());
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

    const bool is_compressible = sim.material_model->is_compressible();

    // Assemble and store the diagonal of the GMG level matrices derived from:
    // 2*eta*(symgrad u, symgrad v) - (if compressible) 2*eta/3*(div u, div v)
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
            FEValues<dim> fe_values_projection (*(sim.mapping),
                                                fe_projection,
                                                quadrature_formula,
                                                update_values);

            const unsigned int   dofs_per_cell   = fe_v.dofs_per_cell;
            const unsigned int   n_q_points      = quadrature_formula.size();

            FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

            std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
            const FEValuesExtractors::Vector velocities (0);

            std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
            std::vector<double> div_phi_u (dofs_per_cell);

            AffineConstraints<double> boundary_constraints;
            boundary_constraints.reinit(locally_relevant_dofs);
            boundary_constraints.add_lines (mg_constrained_dofs_A_block.get_refinement_edge_indices(level));
            boundary_constraints.add_lines (mg_constrained_dofs_A_block.get_boundary_indices(level));
            // let Dirichlet values win over no normal flux:
            boundary_constraints.merge(mg_constrained_dofs_A_block.get_user_constraint_matrix(level),
                                       AffineConstraints<double>::left_object_wins);
            boundary_constraints.close();

            typename DoFHandler<dim>::level_cell_iterator cell = dof_handler_v.begin(level),
                                                          endc = dof_handler_v.end(level);
            for (; cell!=endc; ++cell)
              if (cell->level_subdomain_id() == sim.triangulation.locally_owned_subdomain())
                {
                  cell_matrix = 0;
                  fe_values.reinit (cell);

                  typename DoFHandler<dim>::level_cell_iterator DG_cell(&(sim.triangulation),
                                                                        level,
                                                                        cell->index(),
                                                                        &dof_handler_projection);
                  std::vector<types::global_dof_index> dg_dof_indices(dof_handler_projection.get_fe(0).dofs_per_cell);
                  DG_cell->get_active_or_mg_dof_indices(dg_dof_indices);

                  // For DGQ1, project viscosity from DoF vector to quadrature.
                  std::vector<GMGNumberType> visc_on_quad(n_q_points);
                  if (dof_handler_projection.get_fe().degree == 1)
                    {
                      fe_values_projection.reinit(DG_cell);
                      fe_values_projection.get_function_values(level_viscosity_vector[level],
                                                               dg_dof_indices,
                                                               visc_on_quad);
                    }

                  for (unsigned int q=0; q<n_q_points; ++q)
                    {
                      const double viscosity = (dof_handler_projection.get_fe().degree == 0
                                                ?
                                                level_viscosity_vector[level](dg_dof_indices[0])
                                                :
                                                visc_on_quad[q]);

                      for (unsigned int k=0; k<dofs_per_cell; ++k)
                        {
                          symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);

                          if (is_compressible)
                            div_phi_u[k] = fe_values[velocities].divergence (k, q);
                        }

                      const double JxW = fe_values.JxW(q);
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        for (unsigned int j=0; j<dofs_per_cell; ++j)
                          {
                            cell_matrix(i,j) += 2. * viscosity * (symgrad_phi_u[i]*symgrad_phi_u[j])
                                                * JxW;

                            if (is_compressible)
                              cell_matrix(i,j) +=  (-2./3.) * viscosity * (div_phi_u[i]*div_phi_u[j])
                                                   * JxW;
                          }
                    }

                  cell->get_mg_dof_indices (local_dof_indices);

                  boundary_constraints.distribute_local_to_global (cell_matrix,
                                                                   local_dof_indices,
                                                                   diagonal_matrix);
                }

            diagonal_matrix.compress(VectorOperation::add);

            dealii::LinearAlgebra::distributed::Vector<GMGNumberType> diagonal;
            // This assignment converts from type double to GMGNumberType (float).
            diagonal = diagonal_matrix.get_vector();
            mg_matrices_A_block[level].set_diagonal(diagonal);
          }
        else
          {
            mg_matrices_A_block[level].compute_diagonal();
          }

        // This vector is no longer needed. Resize to 0.
        level_viscosity_vector[level].reinit(0);
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
  const MGTransferMatrixFree<dim,GMGNumberType> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_mg_transfer_A() const
  {
    return mg_transfer_A_block;
  }


  template <int dim, int velocity_degree>
  const MGTransferMatrixFree<dim,GMGNumberType> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_mg_transfer_S() const
  {
    return mg_transfer_Schur_complement;
  }


  template <int dim, int velocity_degree>
  const Table<2, VectorizedArray<double>> &
                                       StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_active_viscosity_table() const
  {
    return active_viscosity_table;
  }


  template <int dim, int velocity_degree>
  const MGLevelObject<Table<2, VectorizedArray<GMGNumberType>>> &
  StokesMatrixFreeHandlerImplementation<dim, velocity_degree>::get_level_viscosity_tables() const
  {
    return level_viscosity_tables;
  }




// explicit instantiation of the functions we implement in this file
#define INSTANTIATE(dim) \
  template class StokesMatrixFreeHandler<dim>; \
  template class StokesMatrixFreeHandlerImplementation<dim,2>; \
  template class StokesMatrixFreeHandlerImplementation<dim,3>;

  ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
}
