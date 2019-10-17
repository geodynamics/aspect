/*
  Copyright (C) 2018 by the authors of the ASPECT code.
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


#ifndef _aspect_mesh_deformation_fast_scape_h
#define _aspect_mesh_deformation_fast_scape_h

#include <aspect/mesh_deformation/interface.h>

/*
 * Define FastScape functions as C functions. Must use exact same function/variable name
 * and type as used in FastScape. All function names must be made lowercase, and an
 * underscore added at the end. Types must be defined as pointers, and sent to
 * fastscape as a reference.
 */

#ifdef __cplusplus
extern"C" {
#endif
//Functions to initialize FastScape
void fastscape_init_();
void fastscape_set_nx_ny_(const int *nnx, const int *nny);
void fastscape_setup_();

//TODO: These were originally double and not const double, make sure this isn't an issue.
void fastscape_set_xl_yl_(const double *xxl,const double *yyl);
void fastscape_set_dt_(double *dtt);
void fastscape_init_h_(double *hp);
void fastscape_set_erosional_parameters_(double *kkf,const double *kkfsed,const double *mm,const double *nnn,
                                         double *kkd,const double *kkdsed,const double *gg1,const double *gg2,const double *pp);

void fastscape_view_();
void fastscape_set_bc_(const int *jbc);
void fastscape_set_v_(double *ux, double *uy);
void fastscape_set_u_(double *up);
void fastscape_set_h_(double *hp);

//Functions to run FastScape
void fastscape_get_step_(int *sstep);
void fastscape_execute_step_();
void fastscape_named_vtk_(double *fp, const double *vexp, int *astep, const char *c, int *length);
void fastscape_copy_h_(double *hp);

//end run
void fastscape_debug_();
void fastscape_destroy_();
#ifdef  __cplusplus
}
#endif

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {
    template<int dim>
    class FastScape : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:

        //FastScape();

        virtual void initialize ();

        virtual
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 ConstraintMatrix &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const;


        /**
         * Declare parameters for the free surface handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm);

      private:
        /**
         * A function object representing the mesh deformation.
         */
        int nstep;
        int bc;
        int array_size;
        double max_timestep;
        double vexp;
        int nx;
        int ny;
        int nz;
        unsigned int additional_refinement;
        unsigned int initial_global_refinement;
        int numx;
        int numy;
        unsigned int bottom;
        unsigned int right;
        unsigned int top;
        unsigned int left;
        double x_extent;
        double y_extent;
        bool slice;
        int fs_seed;
        int surface_res;
        mutable bool restart;
        mutable int restart_step;

        double m;
        double n;
        double p;
        double g;
        double gsed;
        double kff;
        double kfsed;
        double kdd;
        double kdsed;
        double dx;
        double dy;
        double end_time;

        std::array<std::pair<double,double>,dim> grid_extent;
        std::array< unsigned int, dim > table_intervals;
    };
  }
}


#endif
