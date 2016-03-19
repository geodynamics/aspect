/*
 Copyright (C) 2016 by the authors of the ASPECT code.

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
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef __aspect__vofinterface_VOFEngine_h
#define __aspect__vofinterface_VOFEngine_h

// Aspect includes
#include <aspect/global.h>

// Deal II includes
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/numerics/data_out.h>

#include <aspect/vofinterface/QCMid.h>
#include <aspect/vofinterface/VOFUtils.h>

// Definition

namespace aspect
{
  namespace InterfaceTracker
  {
    using namespace dealii;

    template <int dim>
    class VOFEngine
    {
      public:
        VOFEngine ();
        ~VOFEngine ();

        // Component info
        static std::vector<std::string> component_names ();
        static std::vector<DataComponentInterpretation::DataComponentInterpretation>
        component_interpretation ();

        // Code logic
        void init (Function<dim> &func,
                   unsigned int n_samp,
                   double time);
        void calc_normals ();
        void do_step (const LinearAlgebra::BlockVector &new_soln,
                      const LinearAlgebra::BlockVector &old_soln,
                      double timestep);

        // Error estimation
        static std::vector<std::string> error_names ();
        static std::vector<std::string> error_abrev ();
        std::vector<double> calc_error (Function<dim> &func,
                                        unsigned int n_samp,
                                        double time);

        // Setters
        void set_voleps (double new_voleps);
        void set_triangulation (const parallel::distributed::Triangulation<dim> *new_tria);
        void set_mapping (const Mapping<dim> *mapping);
        void set_parent_dofs (const DoFHandler<dim> *new_par_dof);
        void set_comm (const MPI_Comm new_comm_world);
        DoFHandler<dim> *get_dof_handler ();
        Vector<double> &get_state ();
        void clear_dof_handler ();

      private:
        // variable locations
        static const unsigned int n_components = dim + 2;
        static const unsigned int vof_ind = 0;
        static const unsigned int first_normal_ind = 1;
        static const unsigned int d_ind = dim + 1;

        // Configuration vars
        double voleps;

        // Local state vars
        FESystem<dim> interfaceFE;
        DoFHandler<dim> *dof_handler;
        bool normals_calced;
        bool old_vel_set;
        Vector<double> state, deltaState;
        unsigned int dir_order;
        const int qorder;

        // Error vars
        double vol_init;

        // Linking vars
        const parallel::distributed::Triangulation<dim> *triangulation;
        const DoFHandler<dim> *parent_sim_handler;
        const Mapping<dim> *mapping;

        // Communication vars
        MPI_Comm communicator;

        unsigned int world_size;
        unsigned int self_rank;

        // Utility funcs for vol calculations
        void calc_flux (const LinearAlgebra::BlockVector &par_new_soln,
                        const LinearAlgebra::BlockVector &par_old_soln,
                        double timestep,
                        unsigned int dir);
        void update_vof (const LinearAlgebra::BlockVector &par_new_soln,
                         const LinearAlgebra::BlockVector &par_old_soln,
                         double timestep,
                         unsigned int dir);
    };
  }
}

#endif
