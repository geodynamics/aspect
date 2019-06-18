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


#include <aspect/mesh_deformation/fastscape.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/grid/grid_generator.h>


#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {





    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                                    ConstraintMatrix &mesh_velocity_constraints,
                                                                    const std::set<types::boundary_id> &boundary_ids) const
    {

      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      AssertThrow(geometry != nullptr,
                        ExcMessage("Fastscape can only be run with a box model"));
     // AssertThrow(dim == 3,
     //                   ExcMessage("Fastscape can only be applied to 3D runs."));

      if (geometry != nullptr)
        {
          int istep;
          int steps = nstep;
          double xll = geometry->get_extents()[0];
          double yll = geometry->get_extents()[1];

          /*
           * Pointers get allocated only 8 bytes at compile so don't run into
           * allocation issues. This creates an array that can delete itself.
           */
          std::unique_ptr<double[]> h (new double[array_size]);
          std::unique_ptr<double[]> kf (new double[array_size]);
          std::unique_ptr<double[]> kd (new double[array_size]);

          //TODO: Make sure this doesn't run on first timestep.
          double a_dt = 100001; //this->get_timestep();
          double f_dt = a_dt/steps;

          while (f_dt>max_timestep)
          {
        	  steps=steps*2;
        	  f_dt = a_dt/steps;
          }
    	  std::cout<<"Fastscape timestep: "<<f_dt<<std::endl;

          //Currently initializing h this way until we can get info from aspect.
          //TODO: Get h (topography) and later ux/uy from aspect (do we need uplift rate as well
          //for vertical projection?).
          for (int i=0; i<array_size; i++)
            {
              h[i]=0;
              kf[i]=kff;
              kd[i]=kdd;
              if (i>(array_size/2))
                h[i]=1000;
            }

          //Initialize fastscape
          fastscape_init_();
          fastscape_set_nx_ny_(&nx,&ny);
          fastscape_setup_();

          double xl = 10e3;
          double yl = 20e3;
          //set x and y extent
          fastscape_set_xl_yl_(&xl,&yl);

          double dt = 2e3;
          //Set time step
          fastscape_set_dt_(&dt);

          //Initialize topography
          fastscape_init_h_(h.get());

          //Set erosional parameters
          fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);

          //set boundary conditions
          fastscape_set_bc_(&bc);

          //Initialize first time step
          fastscape_get_step_(&istep);

          double summ = 0;
          //loop on time stepping
          //TODO only run this with one process
          do
            {
              //execute step, this increases timestep counter
              fastscape_execute_step_();
              //get value of time step counter
              fastscape_get_step_(&istep);
              //outputs new h values
              fastscape_copy_h_(h.get());
              //output vtk
              fastscape_vtk_(h.get(), &vexp);

              double minnn = h[0];
              double maxxx = h[0];
              summ=0;
              for (int i=0; i<array_size; i++)
                {
                  summ += h[i];
                  if (h[i+1] > maxxx)
                    maxxx = h[i];

                  if (h[i+1] < minnn)
                    minnn = h[i];
                }

              std::cout<<"step: "<<istep<<" h range: "<<minnn<<"  "<<summ/array_size<<"  "<<maxxx<<std::endl;

            }
          while (istep<nstep);

          //output timing
          fastscape_debug_();

          //end FastScape run*/
          fastscape_destroy_();
        }
    }



    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscape");
        {
          prm.declare_entry("X grid dimension", "0",
                            Patterns::Integer(),
                            "Mesh resolution in x direction");
          prm.declare_entry("Y grid dimension", "0",
                            Patterns::Integer(),
                            "Mesh resolution in y direction");
          prm.declare_entry("Number of steps", "0",
                            Patterns::Integer(),
	                        "Number of steps per ASPECT timestep");
		  prm.declare_entry("Maximum timestep", "10e3",
					         Patterns::Double(0),
                            "Maximum timestep for fastscape.");
          prm.declare_entry("Boundary conditions", "1111",
                            Patterns::Integer(),
                            "Boundary conditions where 0 is reflective and 1 is fixed  "
                            "height. Must be given in four digits, where the order is bottom,"
                            "right, top, left.");
          prm.declare_entry("Vertical exaggeration", "3",
                            Patterns::Double(),
                            "Vertical exaggeration for fastscape's VTK file.");

          prm.enter_subsection ("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.5",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Multi-direction slope exponent", "1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Bedrock deposition coefficient", "0",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Sediment deposition coefficient", "0",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Bedrock river incision rate", "2e-6",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Bedrock diffusivity", "1e-1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void FastScape<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          nx = prm.get_integer("X grid dimension");
          ny = prm.get_integer("Y grid dimension");
          nstep = prm.get_integer("Number of steps");
          max_timestep = prm.get_double("Maximum timestep");
          bc = prm.get_integer("Boundary conditions");
          vexp = prm.get_double("Vertical exaggeration");
          array_size = nx*ny;

          prm.enter_subsection("Erosional parameters");
          {
            m = prm.get_double("Drainage area exponent");
            n = prm.get_double("Slope exponent");
            kfsed = prm.get_double("Sediment river incision rate");
            kff = prm.get_double("Bedrock river incision rate");
            kdsed = prm.get_double("Sediment diffusivity");
            kdd = prm.get_double("Bedrock diffusivity");
            g = prm.get_double("Bedrock deposition coefficient");
            gsed = prm.get_double("Sediment deposition coefficient");
            p = prm.get_double("Multi-direction slope exponent");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScape,
                                           "fastscape",
                                           "A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
