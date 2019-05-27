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


#ifndef _aspect_stokes_matrix_free_h
#define _aspect_stokes_matrix_free_h

#include <aspect/global.h>

#include <aspect/simulator.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>

namespace aspect
{
  using namespace dealii;

  /**
   * This namespace contains all matrix-free operators used in the Stokes solver.
   */
  namespace MatrixFreeStokesOperators
  {
    /**
     * Operator for the entire Stokes block.
     */
    template <int dim, int degree_v, typename number>
    class StokesOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number> >
    {
      public:

        /**
         * Constructor.
         */
        StokesOperator ();

        /**
         * Reset the viscosity table.
         */
        void clear ();

        /**
         * Fills in the viscosity table and set the value for the pressure scaling constant.
         */
        void fill_viscosities_and_pressure_scaling(const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                                                   const double pressure_scaling,
                                                   const Triangulation<dim> &tria,
                                                   const DoFHandler<dim> &dof_handler_for_projection);

        /**
         * Returns the viscosity table.
         */
        const Table<2, VectorizedArray<number> > &
        get_viscosity_x_2_table();

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        virtual void compute_diagonal ();

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        virtual void apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                                const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                          const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * Table which stores the viscosity on each quadrature point.
         */
        Table<2, VectorizedArray<number> > viscosity_x_2;

        /**
         * Pressure scaling constant.
         */
        double pressure_scaling;
    };

    template <int dim, int degree_v, typename number>
    StokesOperator<dim,degree_v,number>::StokesOperator ()
      :
      MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number> >()
    {}

    template <int dim, int degree_v, typename number>
    void
    StokesOperator<dim,degree_v,number>::clear ()
    {
      viscosity_x_2.reinit(0, 0);
      MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::BlockVector<number> >::clear();
    }

    template <int dim, int degree_v, typename number>
    void
    StokesOperator<dim,degree_v,number>::
    fill_viscosities_and_pressure_scaling (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                                           const double pressure_scaling,
                                           const Triangulation<dim> &tria,
                                           const DoFHandler<dim> &dof_handler_for_projection)
    {
      FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (*this->data, 0);
      const unsigned int n_cells = this->data->n_macro_cells();
      viscosity_x_2.reinit(n_cells, velocity.n_q_points);

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
            for (unsigned int q=0; q<velocity.n_q_points; ++q)
              viscosity_x_2(cell,q)[i] = 2.0*viscosity_values(local_dof_indices[0]);
          }
      this->pressure_scaling = pressure_scaling;
    }

    template <int dim, int degree_v, typename number>
    const Table<2, VectorizedArray<number> > &
    StokesOperator<dim,degree_v,number>::get_viscosity_x_2_table()
    {
      return viscosity_x_2;
    }

    template <int dim, int degree_v, typename number>
    void
    StokesOperator<dim,degree_v,number>
    ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                   dealii::LinearAlgebra::distributed::BlockVector<number>       &dst,
                   const dealii::LinearAlgebra::distributed::BlockVector<number> &src,
                   const std::pair<unsigned int, unsigned int>           &cell_range) const
    {
      typedef VectorizedArray<number> vector_t;
      FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (data, 0);
      FEEvaluation<dim,degree_v-1,  degree_v+1,1,  number> pressure (data, 1);

      for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
        {
          velocity.reinit (cell);
          velocity.read_dof_values (src.block(0));
          velocity.evaluate (false,true,false);
          pressure.reinit (cell);
          pressure.read_dof_values (src.block(1));
          pressure.evaluate (true,false,false);

          for (unsigned int q=0; q<velocity.n_q_points; ++q)
            {
              SymmetricTensor<2,dim,vector_t> sym_grad_u =
                velocity.get_symmetric_gradient (q);
              vector_t pres = pressure.get_value(q);
              vector_t div = -trace(sym_grad_u);
              pressure.submit_value   (pressure_scaling*div, q);

              sym_grad_u *= viscosity_x_2(cell,q);

              for (unsigned int d=0; d<dim; ++d)
                sym_grad_u[d][d] -= pressure_scaling*pres;

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
    StokesOperator<dim,degree_v,number>
    ::apply_add (dealii::LinearAlgebra::distributed::BlockVector<number> &dst,
                 const dealii::LinearAlgebra::distributed::BlockVector<number> &src) const
    {
      MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::BlockVector<number> >::
      data->cell_loop(&StokesOperator::local_apply, this, dst, src);
    }

    template <int dim, int degree_v, typename number>
    void
    StokesOperator<dim,degree_v,number>
    ::compute_diagonal ()
    {
      // There is currently no need in the code for the diagonal of the entire stokes
      // block. If needed, one could easily construct based on the diagonal of the A
      // block and append zeros to the end for the number of pressure DoFs.
      Assert(false, ExcNotImplemented());
    }

    /**
     * Operator for the pressure mass matrix used in the block preconditioner
     */
    template <int dim, int degree_p, typename number>
    class MassMatrixOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>
    {
      public:

        /**
         * Constructor
         */
        MassMatrixOperator ();

        /**
         * Reset the viscosity table.
         */
        void clear ();

        /**
         * Fills in the viscosity table and set the value for the pressure scaling constant.
         */
        void fill_viscosities_and_pressure_scaling (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                                                    const double pressure_scaling,
                                                    const Triangulation<dim> &tria,
                                                    const DoFHandler<dim> &dof_handler_for_projection);


        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        virtual void compute_diagonal ();

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        virtual void apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                                const dealii::LinearAlgebra::distributed::Vector<number> &src) const;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::Vector<number> &dst,
                          const dealii::LinearAlgebra::distributed::Vector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;


        /**
         * Computes the diagonal contribution from a cell matrix.
         */
        void local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                                     dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                                     const unsigned int                               &dummy,
                                     const std::pair<unsigned int,unsigned int>       &cell_range) const;

        /**
         * Table which stores the viscosity on each quadrature point.
         */
        Table<2, VectorizedArray<number> > one_over_viscosity;

        /**
         * Pressure scaling constant.
         */
        double pressure_scaling;
    };

    template <int dim, int degree_p, typename number>
    MassMatrixOperator<dim,degree_p,number>::MassMatrixOperator ()
      :
      MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number> >()
    {}

    template <int dim, int degree_p, typename number>
    void
    MassMatrixOperator<dim,degree_p,number>::clear ()
    {
      one_over_viscosity.reinit(0, 0);
      MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::clear();
    }

    template <int dim, int degree_p, typename number>
    void
    MassMatrixOperator<dim,degree_p,number>::
    fill_viscosities_and_pressure_scaling (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                                           const double pressure_scaling,
                                           const Triangulation<dim> &tria,
                                           const DoFHandler<dim> &dof_handler_for_projection)
    {
      FEEvaluation<dim,degree_p,degree_p+2,1,number> pressure (*this->data, 0);
      const unsigned int n_cells = this->data->n_macro_cells();
      one_over_viscosity.reinit(n_cells, pressure.n_q_points);

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
            for (unsigned int q=0; q<pressure.n_q_points; ++q)
              one_over_viscosity(cell,q)[i] = 1.0/viscosity_values(local_dof_indices[0]);
          }
      this->pressure_scaling = pressure_scaling;
    }

    template <int dim, int degree_p, typename number>
    void
    MassMatrixOperator<dim,degree_p,number>
    ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                   dealii::LinearAlgebra::distributed::Vector<number>       &dst,
                   const dealii::LinearAlgebra::distributed::Vector<number> &src,
                   const std::pair<unsigned int, unsigned int>           &cell_range) const
    {
      FEEvaluation<dim,degree_p,degree_p+2,1,number> pressure (data);

      for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
        {
          AssertDimension(one_over_viscosity.size(0), data.n_macro_cells());
          AssertDimension(one_over_viscosity.size(1), pressure.n_q_points);

          pressure.reinit (cell);
          pressure.read_dof_values(src);
          pressure.evaluate (true, false);
          for (unsigned int q=0; q<pressure.n_q_points; ++q)
            pressure.submit_value(one_over_viscosity(cell,q)*pressure_scaling*pressure_scaling*
                                  pressure.get_value(q),q);
          pressure.integrate (true, false);
          pressure.distribute_local_to_global (dst);
        }
    }

    template <int dim, int degree_p, typename number>
    void
    MassMatrixOperator<dim,degree_p,number>
    ::apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                 const dealii::LinearAlgebra::distributed::Vector<number> &src) const
    {
      MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::
      data->cell_loop(&MassMatrixOperator::local_apply, this, dst, src);
    }

    template <int dim, int degree_p, typename number>
    void
    MassMatrixOperator<dim,degree_p,number>
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
    MassMatrixOperator<dim,degree_p,number>
    ::local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                              dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                              const unsigned int &,
                              const std::pair<unsigned int,unsigned int>       &cell_range) const
    {
      FEEvaluation<dim,degree_p,degree_p+2,1,number> pressure (data, 0);
      for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
        {
          pressure.reinit (cell);
          AlignedVector<VectorizedArray<number> > diagonal(pressure.dofs_per_cell);
          for (unsigned int i=0; i<pressure.dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<pressure.dofs_per_cell; ++j)
                pressure.begin_dof_values()[j] = VectorizedArray<number>();
              pressure.begin_dof_values()[i] = make_vectorized_array<number> (1.);

              pressure.evaluate (true,false,false);
              for (unsigned int q=0; q<pressure.n_q_points; ++q)
                pressure.submit_value(one_over_viscosity(cell,q)*pressure_scaling*pressure_scaling*
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
     * Operator for the A block of the Stokes matrix. The same class is used for both
     * active and level mesh operators.
     */
    template <int dim, int degree_v, typename number>
    class ABlockOperator
      : public MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number>>
    {
      public:

        /**
         * Constructor
         */
        ABlockOperator ();

        /**
         * Reset the viscosity table.
         */
        void clear ();

        /**
         * Fills in the viscosity table.
         */
        void fill_viscosities(const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                              const Triangulation<dim> &tria,
                              const DoFHandler<dim> &dof_handler_for_projection,
                              const bool for_mg);

        /**
         * Computes the diagonal of the matrix. Since matrix-free operators have not access
         * to matrix elements, we must apply the matrix-free operator to the unit vectors to
         * recover the diagonal.
         */
        virtual void compute_diagonal ();

      private:

        /**
         * Performs the application of the matrix-free operator. This function is called by
         * vmult() functions MatrixFreeOperators::Base.
         */
        virtual void apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                                const dealii::LinearAlgebra::distributed::Vector<number> &src) const;

        /**
         * Defines the application of the cell matrix.
         */
        void local_apply (const dealii::MatrixFree<dim, number> &data,
                          dealii::LinearAlgebra::distributed::Vector<number> &dst,
                          const dealii::LinearAlgebra::distributed::Vector<number> &src,
                          const std::pair<unsigned int, unsigned int> &cell_range) const;

        /**
         * Computes the diagonal contribution from a cell matrix.
         */
        void local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                                     dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                                     const unsigned int                               &dummy,
                                     const std::pair<unsigned int,unsigned int>       &cell_range) const;

        /**
         * Table which stores the viscosity on each quadrature point.
         */
        Table<2, VectorizedArray<number> > viscosity_x_2;

    };

    template <int dim, int degree_v, typename number>
    ABlockOperator<dim,degree_v,number>::ABlockOperator ()
      :
      MatrixFreeOperators::Base<dim, dealii::LinearAlgebra::distributed::Vector<number> >()
    {}

    template <int dim, int degree_v, typename number>
    void
    ABlockOperator<dim,degree_v,number>::clear ()
    {
      viscosity_x_2.reinit(0, 0);
      MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::clear();
    }

    template <int dim, int degree_v, typename number>
    void
    ABlockOperator<dim,degree_v,number>::
    fill_viscosities (const dealii::LinearAlgebra::distributed::Vector<number> &viscosity_values,
                      const Triangulation<dim> &tria,
                      const DoFHandler<dim> &dof_handler_for_projection,
                      const bool for_mg)
    {
      FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (*this->data, 0);
      const unsigned int n_cells = this->data->n_macro_cells();
      viscosity_x_2.reinit(n_cells, velocity.n_q_points);

      std::vector<types::global_dof_index> local_dof_indices(dof_handler_for_projection.get_fe().dofs_per_cell);
      for (unsigned int cell=0; cell<n_cells; ++cell)
        for (unsigned int i=0; i<this->get_matrix_free()->n_components_filled(cell); ++i)
          {

            if (for_mg)
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
            for (unsigned int q=0; q<velocity.n_q_points; ++q)
              viscosity_x_2(cell,q)[i] = 2.0*viscosity_values(local_dof_indices[0]);
          }
    }

    template <int dim, int degree_v, typename number>
    void
    ABlockOperator<dim,degree_v,number>
    ::local_apply (const dealii::MatrixFree<dim, number>                 &data,
                   dealii::LinearAlgebra::distributed::Vector<number>       &dst,
                   const dealii::LinearAlgebra::distributed::Vector<number> &src,
                   const std::pair<unsigned int, unsigned int>           &cell_range) const
    {
      FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (data,0);

      for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
        {
          AssertDimension(viscosity_x_2.size(0), data.n_macro_cells());
          AssertDimension(viscosity_x_2.size(1), velocity.n_q_points);

          velocity.reinit (cell);
          velocity.read_dof_values(src);
          velocity.evaluate (false, true, false);
          for (unsigned int q=0; q<velocity.n_q_points; ++q)
            {
              velocity.submit_symmetric_gradient
              (viscosity_x_2(cell,q)*velocity.get_symmetric_gradient(q),q);
            }
          velocity.integrate (false, true);
          velocity.distribute_local_to_global (dst);
        }
    }

    template <int dim, int degree_v, typename number>
    void
    ABlockOperator<dim,degree_v,number>
    ::apply_add (dealii::LinearAlgebra::distributed::Vector<number> &dst,
                 const dealii::LinearAlgebra::distributed::Vector<number> &src) const
    {
      MatrixFreeOperators::Base<dim,dealii::LinearAlgebra::distributed::Vector<number> >::
      data->cell_loop(&ABlockOperator::local_apply, this, dst, src);
    }

    template <int dim, int degree_v, typename number>
    void
    ABlockOperator<dim,degree_v,number>
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
    ABlockOperator<dim,degree_v,number>
    ::local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                              dealii::LinearAlgebra::distributed::Vector<number>  &dst,
                              const unsigned int &,
                              const std::pair<unsigned int,unsigned int>       &cell_range) const
    {
      FEEvaluation<dim,degree_v,degree_v+1,dim,number> velocity (data, 0);
      for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
        {
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
                  velocity.submit_symmetric_gradient
                  (viscosity_x_2(cell,q)*velocity.get_symmetric_gradient(q),q);
                }
              velocity.integrate (false,true);

              diagonal[i] = velocity.begin_dof_values()[i];
            }

          for (unsigned int i=0; i<velocity.dofs_per_cell; ++i)
            velocity.begin_dof_values()[i] = diagonal[i];
          velocity.distribute_local_to_global (dst);
        }
    }
  }



  template<int dim>
  class StokesMatrixFreeHandler
  {
    public:
      /**
       * Initialize this class, allowing it to read in
       * relevant parameters as well as giving it a reference to the
       * Simulator that owns it, since it needs to make fairly extensive
       * changes to the internals of the simulator.
       */
      StokesMatrixFreeHandler(Simulator<dim> &, ParameterHandler &prm);

      /**
       * Destructor.
       */
      ~StokesMatrixFreeHandler();

      /**
       * Solves the Stokes linear system matrix-free. This is called
       * by Simulator<dim>::solve_stokes().
       */
      std::pair<double,double> solve();

      /**
       * Allocates and sets up the members of the StokesMatrixFreeHandler. This
       * is called by Simulator<dim>::setup_dofs()
       */
      void setup_dofs();

      /**
       * Evalute the MaterialModel to query for the viscosity on the active cells,
       * project this viscosity to the multigrid hierarchy, and cache the information
       * for later usage.
       */
      void evaluate_viscosity();

      /**
       * Get the workload imbalance of the distribution
       * of the level hierarchy.
       */
      double get_workload_imbalance();

      /**
       * Add correction to system RHS for non-zero boundary condition.
       */
      void correct_stokes_rhs();

      /**
       * Declare parameters. (No actual parameters at the moment).
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Parse parameters. (No actual parameters at the moment).
       */
      void parse_parameters (ParameterHandler &prm);


    private:

      Simulator<dim> &sim;

      DoFHandler<dim> dof_handler_v;
      DoFHandler<dim> dof_handler_p;
      DoFHandler<dim> dof_handler_projection;

      FESystem<dim> stokes_fe;
      FESystem<dim> fe_v;
      FESystem<dim> fe_p;
      FESystem<dim> fe_projection;

      // TODO: velocity degree not only 2, Choosing quadrature degree?
      typedef MatrixFreeStokesOperators::StokesOperator<dim,2,double> StokesMatrixType;
      typedef MatrixFreeStokesOperators::MassMatrixOperator<dim,1,double> MassMatrixType;
      typedef MatrixFreeStokesOperators::ABlockOperator<dim,2,double> ABlockMatrixType;

      StokesMatrixType stokes_matrix;
      ABlockMatrixType velocity_matrix;
      MassMatrixType mass_matrix;

      ConstraintMatrix constraints_v;
      ConstraintMatrix constraints_p;
      ConstraintMatrix constraints_projection;

      MGLevelObject<ABlockMatrixType> mg_matrices;
      MGConstrainedDoFs              mg_constrained_dofs;
      MGConstrainedDoFs mg_constrained_dofs_projection;

      dealii::LinearAlgebra::distributed::Vector<double> active_coef_dof_vec;
      MGLevelObject<dealii::LinearAlgebra::distributed::Vector<double> > level_coef_dof_vec;


      MGTransferMatrixFree<dim,double> mg_transfer;
  };
}


#endif
