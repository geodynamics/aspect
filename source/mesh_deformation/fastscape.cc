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
#include <aspect/utilities.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>

#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    void
    FastScape<dim>::initialize ()
    {
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      AssertThrow(geometry != nullptr,
                  ExcMessage("Fastscape can only be run with a box model"));
      AssertThrow(dim == 3,
                  ExcMessage("Fastscape can only be applied to 3D runs."));

      x_extent = geometry->get_extents()[0];
      y_extent = geometry->get_extents()[1];

      //second is for maximum coordiantes, first for minimum.
      grid_extent[0].second = geometry->get_extents()[0];
      grid_extent[1].second = geometry->get_extents()[1];
      grid_extent[0].first = 0;
      grid_extent[1].first = 0;


      //TODO: There has to be a better type to use to get this.
      const std::pair<int, int > repetitions = geometry->get_repetitions();
      int x_repetitions = repetitions.first;
      int y_repetitions = repetitions.second;

      /*
       * Set nx and ny based on repetitions and chosen max refinement.
       * TODO: There needs to be a good error message if refinement is chosen
       * incorrectly. Right now it just throws a segmentation fault.
       */
      nx = 1+std::pow(2,refinement)*x_repetitions;
      ny = 1+std::pow(2,refinement)*y_repetitions;
      array_size = nx*ny-1;
      table_intervals[0] = nx;
      table_intervals[1] = ny;
    }

    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             ConstraintMatrix &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_ids) const
    {
    	if (this->get_timestep_number()!=0)
    	{
      int istep;
      int steps = nstep;

      /*
       * Pointers get allocated only 8 bytes at compile so don't run into
       * allocation issues. This creates an array that can delete itself.
       */
      std::unique_ptr<double[]> h_temp (new double[array_size]);
      std::unique_ptr<double[]> posx_temp (new double[array_size]);
      std::unique_ptr<double[]> posy_temp (new double[array_size]);

      std::unique_ptr<double[]> vx_temp (new double[array_size]);
      std::unique_ptr<double[]> vy_temp (new double[array_size]);
      std::unique_ptr<double[]> vz_temp (new double[array_size]);
      std::unique_ptr<double[]> V (new double[array_size]);


      //Set h_temp to zero everywhere.
      for (int i=0; i<=array_size; i++)
      {
        h_temp[i]=0;
        vx_temp[i]=0;
        vy_temp[i]=0;
        vz_temp[i]=0;
        V[i] = 0;
      }

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

      // Get a quadrature rule that exists only on the corners
      QTrapez<dim-1> face_corners;

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        face_corners,
                                        update_values |
                                        update_quadrature_points);

      // loop over all of the surface cells and if one less than h/3 away from
      // the top or bottom surface, evaluate the density on that face
      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();

       for (; cell != endc; ++cell)
         if (cell->is_locally_owned() && cell->at_boundary())
           for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
             if (cell->face(face_no)->at_boundary())
               {
                 if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                   continue;

                 std::vector<Tensor<1,dim> > vel( face_corners.size() );
                 fe_face_values.reinit( cell, face_no);
                 fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel );

                 //Find out dx and dy, and the total number of x points per row.
                 double dx = cell->extent_in_direction(0);
                 double dy = cell->extent_in_direction(1);
                 double numx = 1+x_extent/dx;

                 for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                   {
                     const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                     //const double vertexx = vel[corner][dim-1];
                     //const double elevation = vertex(dim-1);

                     double indx = 1+vertex(0)/dx;
                     double indy = 1+vertex(1)/dy;
                     double index = (indy-1)*numx+indx;

                     posx_temp[index-1] = vertex(0);
                     posy_temp[index-1] = vertex(1);
                     h_temp[index-1] = vertex(dim-1);

                     vx_temp[index-1] = vel[corner][0];
                     vy_temp[index-1] = vel[corner][1];
                     vz_temp[index-1] = vel[corner][2];

                   }
               }

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          std::unique_ptr<double[]> h (new double[array_size]);
          std::unique_ptr<double[]> posx (new double[array_size]);
          std::unique_ptr<double[]> posy (new double[array_size]);

          std::unique_ptr<double[]> vx (new double[array_size]);
          std::unique_ptr<double[]> vy (new double[array_size]);
          std::unique_ptr<double[]> vz (new double[array_size]);
          std::unique_ptr<double[]> kf (new double[array_size]);
          std::unique_ptr<double[]> kd (new double[array_size]);
          /* Initialize kf and kd across array, and set h values to what proc zero has.
           * TODO: Find a cleaner way to initialize these to the same value,
           * and find how to initialize if they're spatially dependent. Also,
           * is it faster to do this here or earlier?
           */
          for (int i=0; i<=array_size; i++)
            {
              h[i]=h_temp[i];
              posx[i]=posx_temp[i];
              posy[i]=posy_temp[i];

              vx[i]=vx_temp[i];
              vy[i]=vy_temp[i];
              vz[i]=vz_temp[i];
              kf[i]=kff;
              kd[i]=kdd;
            }

          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
            {
              MPI_Status status;
              MPI_Recv(h_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);
              MPI_Recv(posx_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);
              MPI_Recv(posy_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

              MPI_Recv(vx_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);
              MPI_Recv(vy_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);
              MPI_Recv(vz_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

              /*
               * To keep indexing correct, each proc fills h_temp with zeroes in locations that
               * it doesn't own. This checks whether or not we have an actual value.
               */
              for (int i=0; i<=array_size; i++)
                if (h_temp[i]>0)
                {
                  h[i]=h_temp[i];
                  posx[i]=posx_temp[i];
                  posy[i]=posy_temp[i];

                  vx[i]=vx_temp[i];
                  vy[i]=vy_temp[i];
                  vz[i]=vz_temp[i];
                }

            }

          //Set h_temp as initial h values now.
          for (int i=0; i<=array_size; i++)
            {
              h_temp[i] = h[i];
            }

          double a_dt = this->get_old_timestep()/year_in_seconds;

          double f_dt = a_dt/steps;
          while (f_dt>max_timestep)
            {
              steps=steps*2;
              f_dt = a_dt/steps;
            }
          std::cout<<"Fastscape timestep: "<<f_dt<<"  "<<a_dt<<"  "<<steps<<std::endl;

          //Initialize fastscape
          fastscape_init_();
          fastscape_set_nx_ny_(&nx,&ny);
          fastscape_setup_();

          //set x and y extent
          fastscape_set_xl_yl_(&x_extent,&y_extent);

          //Set time step
          f_dt = 1e6;
          fastscape_set_dt_(&f_dt);

          //Initialize topography
          fastscape_init_h_(h.get());

          //Set erosional parameters
          fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);

          //set boundary conditions
          fastscape_set_bc_(&bc);

          //Initialize first time step
          fastscape_get_step_(&istep);
          //fastscape_vtk_(h.get(), &vexp);
          steps =5;

           if(a_dt > 0)
           {
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
               //std::cout<<steps<<"  "<<istep<<std::endl;
             }
           while (istep<steps);
           }

          //output timing
          fastscape_debug_();

          //end FastScape run
          fastscape_destroy_();

          if(a_dt>0)
          {
          for (int i=0; i<=array_size; i++)
        	  V[i] = (h[i] - h_temp[i])/a_dt;
          }


        }
      else
        {
          MPI_Send(h_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
          MPI_Send(posx_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
          MPI_Send(posy_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          MPI_Send(vx_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
          MPI_Send(vy_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
          MPI_Send(vz_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
        }

      int components = 3;
      TableIndices<dim-1> size_idx;
      for (unsigned int d=0; d<dim-1; ++d)
        size_idx[d] = table_intervals[d];

       Table<dim-1,double> data_table;
       data_table.TableBase<dim-1,double>::reinit(size_idx);
       std::vector<Table<dim-1,double> > data_tables(components,data_table);

       TableIndices<dim-1> idx;
       //Indexes through y and then x
       for (unsigned int k=0; k<components; ++k)
       {
    	   std::cout<<k<<std::endl;
       for (unsigned int i=0; i<data_table.size()[1]; ++i)
         {
           idx[1] = i;
           for (unsigned int j=0; j<data_table.size()[0]; ++j)
             {
               idx[0] = j;
               if(k==2)
            	   data_tables[k](idx) = V[nx*i+j];
               else
            	 data_tables[k](idx) = 0;

               std::cout<<data_tables[k](idx)<<"  ";
             }
           std::cout<<"  "<<std::endl;
         }
    }

       /*for (unsigned int i = 0; i < components; i++)
         {
           velocities[i]
             = std_cxx14::make_unique<Functions::InterpolatedUniformGridData<dim-1>> (grid_extent,
                                                                                  table_intervals,
                                                                                  data_tables[i]);
         }*/
    	}
    }


    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscape");
        {
          prm.declare_entry("Number of steps", "5",
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
          prm.declare_entry("Maximum aspect refinement", "1",
                            Patterns::Integer(),
                            "Refinement level expected at surface to determine"
                            "proper nx and ny values");

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
          nstep = prm.get_integer("Number of steps");
          max_timestep = prm.get_double("Maximum timestep");
          bc = prm.get_integer("Boundary conditions");
          vexp = prm.get_double("Vertical exaggeration");
          nstep = prm.get_integer("Maximum aspect refinement");
          refinement = prm.get_integer("Maximum aspect refinement");

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
