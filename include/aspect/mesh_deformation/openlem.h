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


#ifndef _aspect_mesh_deformation_openlem_h
#define _aspect_mesh_deformation_openlem_h

#include <aspect/mesh_deformation/interface.h>
#include <aspect/simulator_access.h>
#include <openlem.cpp>
#include <deal.II/base/parsed_function.h>

#include <aspect/mesh_deformation/openlem.h>

namespace openlem
{
  class Connector
  {
    public:
      double  hscale, x0, y0, alpha, dx0, dy0, dalpha;
      vector<vector<double> >  x, y, vx, vy, vz;
      Grid<>  *g;

      Connector ( Grid<> *g, double hscale = 1, double x0 = 0, double y0 = 0, double alpha = 0 )
      {
        std::cout << "connector construct" << std::endl;
        this->g = g;
        this->hscale = hscale;
        this->x0 = x0;
        this->y0 = y0;
        this->alpha = alpha;
        dx0 = dy0 = dalpha = 0;
        x.resize(g->m);
        for ( int i = 0; i < g->m; ++i )  x[i].resize(g->n);
        y.resize(g->m);
        for ( int i = 0; i < g->m; ++i )  y[i].resize(g->n);
        vx.resize(g->m);
        for ( int i = 0; i < g->m; ++i )  vx[i].resize(g->n);
        vy.resize(g->m);
        for ( int i = 0; i < g->m; ++i )  vy[i].resize(g->n);
        vz.resize(g->m);
        for ( int i = 0; i < g->m; ++i )  vz[i].resize(g->n);

        updateXY();
      }

      void updateXY()
      {
        std::cout << "connector updateXY" << std::endl;
        double  c = cos(alpha);
        double  s = sin(alpha);
        for ( int i = 0; i < g->m; ++i )
          for ( int j = 0; j < g->n; ++j )
            {
              x[i][j] = x0 + hscale*(c*i-s*j);
              y[i][j] = y0 + hscale*(s*i+c*j);
            }
      }

      void convertVelocities()
      {
        std::cout << "connector convertVelocities" << std::endl;
        int     n = 0;
        double  xmean = 0, ymean = 0, vxmean = 0, vymean = 0, cross = 0, rsq = 0;
        double  c = cos(alpha);
        double  s = sin(alpha);

        for ( int i = 0; i < g->m; ++i )
          for ( int j = 0; j < g->n; ++j )
            {
// Convert velocities to OpenLEM coordinate system
              double  tmp = vx[i][j];
              vx[i][j] = (c*vx[i][j]+s*vy[i][j])/hscale;;
              vy[i][j] = (-s*tmp+c*vy[i][j])/hscale;
// Only continental area
              if ( g->getNode(i,j)->b == 0 )
                {
                  ++n;
                  xmean += i;
                  ymean += j;
                  vxmean += vx[i][j];
                  vymean += vy[i][j];
                  cross += i*vy[i][j]-j*vx[i][j];
                  rsq += i*i+j*j;
                }
            }
        assert(n > 0);
        xmean /= n;
        ymean /= n;
        vxmean /= n;
        vymean /= n;
        cross /= n;
        rsq /= n;
        printf ( "%e %e\n", xmean, ymean );
        dalpha = (cross-xmean*vymean+ymean*vxmean)/(rsq-xmean*xmean-ymean*ymean);
        std::cout << "cross = " << cross << ", rsq = " << rsq << ", cross-rsq = " << cross-rsq << ", mean diff top = " << -xmean *vymean+ymean *vxmean << ", mean diff bottom = " << -xmean *xmean-ymean *ymean << std::endl;
        dx0 = vxmean+dalpha*ymean;
        dy0 = vymean-dalpha*xmean;
        printf ( "%e %e %e\n", dalpha, dx0, dy0 );
        //std::cout << " vx:vy = " << 0 << ":" << 0 << " = "<< vx[0][0] << " : " << vy[0][0] << std::endl;
        for ( int i = 0; i < g->m; i=i+10 )
          for ( int j = 0; j < g->n; j=j+10 )
            {
              //if (vx[i][j] > 0 && vy[i][j] > 0)
              std::cout << " vx:vy = " << i << ":" << j << " = "<< vx[i][j] << " : " << vy[i][j] << std::endl;
            }
        for ( int i = 0; i < g->m; ++i )
          for ( int j = 0; j < g->n; ++j )
            {
              vx[i][j] -= dx0-dalpha*j;
              vy[i][j] -= dy0+dalpha*i;
            }

        double dx0_old = dx0;
        dx0 = (c*dx0-s*dy0)*hscale;
        dy0 = (s*dx0_old+c*dy0)*hscale;

        //std::cout << "residual x:y = " << 0 << ":" << 0 << " = "<< vx[0][0] i                     << " : " << vy[0][0] << std::endl;
        //std::cout << "residual x:y = " << 0 << ":" << openlem_ny << " = "<< vx[0][openlem_ny-1]            << " : " << vy[0][openlem_ny-1] << std::endl;
        //std::cout << "residual x:y = " << 0 << ":" << 0 << " = "<< vx[openlem_nx-1][openlem_ny-1] << " : " << vy[openlem_nx-1][openlem_ny-1] << std::endl;
        //std::cout << "residual x:y = " << 0 << ":" << 0 << " = "<< vx[openlem_nx-1][0]            << " : " << vy[openlem_nx-1][0] << std::endl;
        for ( int i = 0; i < g->m; i=i+10 )
          for ( int j = 0; j < g->n; j=j+10 )
            {
              //if (vx[i][j] > 0 && vy[i][j] > 0)
              std::cout << "residual x:y = " << i << ":" << j << " = "<< vx[i][j] << " : " << vy[i][j] << std::endl;
            }
      }

      void computeUpliftRate()
      {
        std::cout << "connector computeUpliftRate" << std::endl;
        for ( int i = 0; i < g->m; ++i )
          for ( int j = 0; j < g->n; ++j )
            if ( !g->getNode(i,j)->b&1 )
              {
                double  dhdx = vx[i][j]<0 ?
                               g->getNodeP(i+1,j)->h-g->getNode(i,j)->h :
                               g->getNode(i,j)->h-g->getNodeP(i-1,j)->h;
                double  dhdy = vy[i][j]<0 ?
                               g->getNodeP(i,j+1)->h-g->getNode(i,j)->h :
                               g->getNode(i,j)->h-g->getNodeP(i,j-1)->h;
                std::cout << "before dhdx:y = " << dhdx << ':' << dhdy << ", u = " << g->getNode(i,j)->u << ", vx:y:z = " << vx[i][j] << ":" << vy[i][j] << ":" << vz[i][j] << std::endl;
                g->getNode(i,j)->u = vz[i][j]-vx[i][j]*dhdx-vy[i][j]*dhdy;
                std::cout << "dhdx:y = " << dhdx << ':' << dhdy << ", u = " << g->getNode(i,j)->u << ", vx:y:z = " << vx[i][j] << ":" << vy[i][j] << ":" << vz[i][j] << std::endl;
              }
      }

      void updateCoordinateSystem ( double dt )
      {
        std::cout << "connector updateCoordinatesystem x0:y0 = " << x0 << ":" << y0 << ", dx0:dy0 = " << dx0 << ":" << dy0 << ", dt = " << dt << ", aplha = " << alpha << ", dalhpa = " << dalpha <<std::endl;
        x0 += dx0*dt;
        y0 += dy0*dt;
        alpha += dalpha*dt;
        std::cout << "connector updateCoordinatesystem x0:y0 = " << x0 << ":" << y0 << ", dx0:dy0 = " << dx0 << ":" << dy0 << ", dt = " << dt << ", alpha = " << alpha << std::endl;
        updateXY();
      }

      double interpolation(double x, double y, const std::vector<std::vector<double>> &interperolat) const
      {

        double c = cos(alpha);
        double s = sin(alpha);
        double x_temp = (x - x0)/hscale;
        double y_temp = (y- y0)/hscale;

        x = c*x_temp+s*y_temp;
        y =-s*x_temp + c*y_temp;
        //std::cout << "new x:y = " << x << ":" << y << std::endl;

        int i = floor(x);
        if (i < 0)
          {
            i = 0;
          }
        else if (i >= interperolat.size())
          {
            i = interperolat.size()-1;
          }
        int j = floor(y);
        if (j < 0)
          {
            j = 0;
          }
        else if (j >= interperolat[i].size())
          {
            j = interperolat[i].size()-1;
          }
        return g->getNode(i,j)->h;//interperolat[i][j];
      }

      void write_vtk(double reference_surface_height, int timestep)
      {
        std::cout << "writing openlem VTK file" << std::endl;
        std::vector<double> grid_x(0);
        std::vector<double> grid_y(0);
        std::vector<double> grid_z(0);
        std::vector<double> grid_elevation(0);

        std::vector<std::vector<size_t> > grid_connectivity(0);

        unsigned int n_cell = x.size() * y.size();
        unsigned int n_p = (x.size() + 1) * (y.size() + 1);


        grid_x.resize(n_p);
        grid_z.resize(n_p);
        grid_y.resize(n_p);

        grid_elevation.resize(n_p);
        grid_connectivity.resize(n_cell,std::vector<size_t>((2-1)*4));

        size_t counter = 0;
        for (size_t j = 0; j <= y.size()-1; ++j)
          {
            for (size_t i = 0; i <= x.size()-1; ++i)
              {
                grid_x[counter] = x[i][j];
                grid_y[counter] = y[i][j];
                grid_z[counter] = reference_surface_height + g->getNode(i,j)->h;
                grid_elevation[counter] = g->getNode(i,j)->h;
                counter++;
              }
          }

        counter = 0;
        for (size_t j = 1; j <= y.size()-1; ++j)
          {
            for (size_t i = 1; i <= x.size()-1; ++i)
              {
                grid_connectivity[counter][0] = i + (j - 1) * (x.size()-1 + 1) - 1;
                grid_connectivity[counter][1] = i + 1 + (j - 1) * (x.size()-1 + 1) - 1;
                grid_connectivity[counter][2] = i + 1  + j * (x.size()-1 + 1) - 1;
                grid_connectivity[counter][3] = i + j * (x.size()-1 + 1) - 1;
                counter++;
              }
          }
        std::stringstream buffer;
        std::ofstream myfile;
        myfile.open ("openlem_surface_" + std::to_string(timestep) +".vtu");
        buffer << "<?xml version=\"1.0\" ?> " << std::endl;
        buffer << R"(<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">)" << std::endl;
        buffer << "<UnstructuredGrid>" << std::endl;
        buffer << "<FieldData>" << std::endl;
        buffer << R"(<DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">0</DataArray>)" << std::endl;
        buffer << "</FieldData>" << std::endl;
        buffer << "<Piece NumberOfPoints=\""<< n_p << "\" NumberOfCells=\"" << n_cell << "\">" << std::endl;
        buffer << "  <Points>" << std::endl;
        buffer << R"(    <DataArray type="Float32" NumberOfComponents="3" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer << grid_x[i] << " " << grid_y[i] << " " << grid_z[i] << std::endl;
          }
        buffer << "    </DataArray>" << std::endl;
        buffer << "  </Points>" << std::endl;
        buffer << std::endl;
        buffer << "  <Cells>" << std::endl;
        buffer << R"(    <DataArray type="Int32" Name="connectivity" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_cell; ++i)
          buffer << grid_connectivity[i][0] << " " <<grid_connectivity[i][1] << " " << grid_connectivity[i][2] << " " << grid_connectivity[i][3] << std::endl;
        buffer << "    </DataArray>" << std::endl;
        buffer << R"(    <DataArray type="Int32" Name="offsets" format="ascii">)" << std::endl;
        for (size_t i = 1; i <= n_cell; ++i)
          buffer << i * 4 << " ";
        buffer << std::endl << "    </DataArray>" << std::endl;
        buffer << R"(    <DataArray type="UInt8" Name="types" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_cell; ++i)
          buffer << "9" << " ";
        buffer <<  std::endl <<"    </DataArray>" << std::endl;
        buffer << "  </Cells>" << std::endl;

        buffer << "  <PointData Scalars=\"scalars\">" << std::endl;

        buffer << R"(<DataArray type="Float32" Name="Elevation" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer <<  grid_elevation[i] << std::endl;
          }
        buffer << "</DataArray>" << std::endl;
        buffer << "  </PointData>" << std::endl;
        //buffer << "  <PointData Vector=\"velocity\">" << std::endl;
        //buffer << R"(<DataArray type="Float32" NumberOfComponents="3" Name="Velocity" format="ascii">)" << std::endl;

        //for (size_t i = 0; i < n_p; ++i)
        //  {
        //    buffer <<  1 << " " << 2 << " " << 3 << std::endl;
        //  }
        //buffer << "</DataArray>" << std::endl;
        buffer << " </Piece>" << std::endl;
        buffer << " </UnstructuredGrid>" << std::endl;
        buffer << "</VTKFile>" << std::endl;
        myfile << buffer.str();
        buffer.str(std::string());

        //std::cout << "                                                                                \r";
        //std::cout.flush();
        std::cout << "Finished writing openlem VTK file" << std::endl;
      }
  };
}
namespace aspect
{
  namespace MeshDeformation
  {
    /**
     * A class that represents a mesh deformation function that can be
     * prescribed on the boundary of the domain.
     *
     * @ingroup MeshDeformation
     */
    template <int dim>
    class OpenLEM : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        struct PositionVelocity
        {
          int x;
          int y;
          Tensor<1,dim> velocity;
        };
        /**
         * Constructor.
         */
        OpenLEM();

        /**
         * Initialize variables for openlem.
         */
        virtual void initialize () override;

        /**
         *
         */
        void update() override;

        /**
         * A function that creates constraints for the velocity of certain mesh
         * vertices (e.g. the surface vertices) for a specific boundary.
         * The calling class will respect
         * these constraints when computing the new vertex positions.
         */
        void
        compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                 AffineConstraints<double> &mesh_velocity_constraints,
                                                 const std::set<types::boundary_id> &boundary_id) const override;

        /**
         * Returns whether or not the plugin requires surface stabilization
         */
        bool needs_surface_stabilization () const override;

        /**
         * Declare parameters for the free surface handling.
         */
        static
        void declare_parameters (ParameterHandler &prm);

        /**
         * Parse parameters for the free surface handling.
         */
        void parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Execute openlem
         */
        void execute_openlem(openlem::Grid<openlem::Node> &grid,
                             //std::vector<double> &elevation,
                             //std::vector<double> &extra_vtk_field,
                             //std::vector<double> &velocity_x,
                             //std::vector<double> &velocity_y,
                             //std::vector<double> &velocity_z,
                             const double &openlem_timestep_in_years,
                             const unsigned int &openlem_iterations) const;

        /**
         * Function to fill the openlem arrays (height and velocities) with the data received from ASPECT in the correct index order.
         */
        void fill_openlem_arrays(openlem::Grid<openlem::Node> &grid,
                                 // std::vector<double> &elevation,
                                 // std::vector<double> &bedrock_transport_coefficient_array,
                                 // std::vector<double> &bedrock_river_incision_rate_array,
                                 // std::vector<double> &velocity_x,
                                 // std::vector<double> &velocity_y,
                                 // std::vector<double> &velocity_z,
                                 std::vector<std::vector<double>> &temporary_variables);

        /**
         * Function to get the ASPECT topography and velocities at the surface, and an index for transferring these to openlem.
         */
        std::vector<std::vector<double>> get_aspect_values() const;

        /**
               * define a old and new grid, so that we can compute a difference
               */
        openlem::Grid<> grid_old;
        openlem::Grid<> grid_new;
        openlem::Connector connector;
        unsigned int openlem_nx;
        unsigned int openlem_ny;
        unsigned int openlem_dx;
        unsigned int openlem_dy;
        unsigned int openlem_x_extent;
        unsigned int openlem_y_extent;
        unsigned int aspect_nx;
        unsigned int aspect_ny;
        unsigned int aspect_dx;
        unsigned int aspect_dy;
        unsigned int aspect_x_extent;
        unsigned int aspect_y_extent;
        std::vector<std::vector<double>> mesh_velocity_z;
        std::vector<std::vector<double>> aspect_mesh_velocity_z;
        /**
         * Variable to hold ASPECT domain extents.
         */
        std::array<std::pair<double,double>,dim> grid_extent;

        /**
         * Table for interpolating openlem surface velocities back to ASPECT.
         */
        std::array<unsigned int, dim> table_intervals;

        /**
         * Whether or not to use the ghost nodes.
         */
        bool use_ghost_nodes;

        /**
         * Check whether openlem needs to be restarted. This is used as
         * a mutable bool because we determine whether the model is being resumed in
         * initialize(), and then after reinitializing openlem we change it to false
         * so it does not initialize openlem again in future timesteps.
         * TODO: There is probably a better way to do this, and restarts should be rolled into
         * the general ASPECT restart.
         */
        mutable bool restart;

        /**
         * How many levels openlem should be refined above the maximum ASPECT surface resolution.
         */
        unsigned int additional_refinement_levels;

        /**
         * Maximum expected refinement level at ASPECT's surface.
         * This and resolution_difference are required to properly transfer node data from
         * ASPECT to openlem.
         */
        unsigned int maximum_surface_refinement_level;

        /**
         * User set openlem Y extent for a 2D ASPECT model.
         */
        double openlem_y_extent_2d;

        /**
         * A time (in seconds) at which the last graphical output was supposed
         * to be produced. Used to check for the next necessary output time.
         */
        mutable double last_output_time;

        /**
         * Suggestion for the number of openlem steps to run for every ASPECT timestep,
         * where the openlem timestep is determined by ASPECT_timestep_length divided by
         * this parameter.
         */
        unsigned int openlem_steps_per_aspect_step;

        /**
         * Maximum timestep allowed for openlem, if the suggested timestep exceeds this
         * limit it is repeatedly divided by 2 until the final timestep is smaller than this parameter.
         */
        double maximum_openlem_timestep;

        /**
         * Difference in refinement levels expected at the ASPECT surface,
         * where this would be set to 2 if 3 refinement levels are set at the surface.
         * This and surface_resolution are required to properly transfer node data from
         * ASPECT to FastScape.
         *
         * TODO: Should this be kept this way, or make it so the input is the expected levels
         * of refinement at the surface, and we can subtract one within the code? Also,
         * it would be good to find a way to check these are correct, because they are a
         * common source of errors.
         */
        unsigned int surface_refinement_difference;


        /**
         * Node tolerance for how close a ASPECT node must be to the FastScape node
         * for the value to be transferred. This is only necessary if use_v is set to 0
         * and the free surface is used to advect the surface with a normal projection, or
         * if there is a surface refinement level difference leading to excess interpolation
         * points in areas of high ASPECT resolution.
         */
        double node_tolerance;

        /**
        * The velocity interpolation data
        */
        Table<dim,double> velocity_table;
    };
  }
}


#endif
