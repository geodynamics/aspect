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
#include <deal.II/base/exceptions.h>
#include <limits>
#include <openlem.cpp>
#include <deal.II/base/parsed_function.h>

//#include <aspect/mesh_deformation/openlem.h>
#include <string>

namespace openlem
{
  class OceanGrid : public Grid<Node>
  {
    public:
      double  odiff;

      OceanGrid ( int m = 1, int n = 1 ) : Grid<Node>(m,n)
      {
        addKey ( "od", "od", &odiff, sizeof(odiff), 8, 0 );
      }

      void clearOceans()
      {
        for ( int i = 0; i < m; ++i )
          for ( int j = 0; j < n; ++j )  getNode(i,j)->b &= 1;
      }

      void markOcean ( Point p, double l )
      {
        if ( getNode(p)->h < l )
          {
            getNode(p)->l = l;
            getNode(p)->b |= 2;
            vector<Point>  neigh = getNeighbors(p);
            for ( vector<Point>::iterator u = neigh.begin(); u != neigh.end(); ++u )
              if ( getNode(u)->drainsTo(p) )
                markOcean(*u,l);
          }
      }

      vector<PointValue<double>> findDeltas ( double dt )
      {
        vector <PointValue<double>>  delta;
        for ( int i = 0; i < m; ++i )
          for ( int j = 0; j < n; ++j )
            if ( getNode(i,j)->b & 2 )
              {
                if ( getNode(i,j)->qs )
                  delta.push_back(PointValue<double>(Point(i,j),getNode(i,j)->qs));
                getNode(i,j)->qp = getNode(i,j)->h+getNode(i,j)->u*dt;
              }
        sort<double>(delta);
        for ( int i = 0; i < delta.size(); ++i )
          printf ( "%i %i %f %e\n", delta[i].p.i, delta[i].p.j, getNode(delta[i].p)->h, delta[i].d );
        if ( offset.size() < m*n )
          {
            printf ( "Computing offsets\n" );
            offset.resize(m*n);
            int  k = 0;
            for ( int i = 0; i < m; ++i )
              for ( int j = 0; j < n; ++j )
                {
                  offset[k].p = Point(i,j);
                  int di = i <= m/2 ? i : m-i;
                  int dj = j <= n/2 ? j : n-j;
                  offset[k++].d = di*di+dj*dj;
                }
            offset[0].d = -1;
            sort(offset);
          }
// Use qp for sediment surface
        for ( int i = 0; i < delta.size(); ++i )
          {
            double v = delta[i].d*dt;
            double tmp;
            double vd = 0;
            Point p = delta[i].p;
            vector<PointValue<float>>::iterator  it = offset.begin();
            while ( it != offset.end() )
              {
                if ( (tmp=getNode(p)->qp+v) <= getNode(p)->l )
                  {
                    getNode(p)->qp = tmp;
                    vd += getNode(p)->qp-getNode(p)->h;
                    break;
                  }
                else
                  {
                    v = v-(getNode(p)->l-getNode(p)->qp);
                    getNode(p)->qp = getNode(p)->l;
                    do
                      {
                        p = (*(it++)).p;
                        p = Point((delta[i].p.i+p.i)%m,(delta[i].p.j+p.j)%n);
                      }
                    while ( (getNode(p)->b&2) == 0 );
                  }
              }
          }
        for ( int i = 0; i < m; ++i )
          for ( int j = 0; j < n; ++j )
            {
              Node  *pn = getNode(i,j);
//        pn->l = 0;
              if ( pn->b&2 )
                {
                  pn->qp = (pn->qp-pn->h)/dt;
//          pn->l = pn->qp;
                }
            }


        return delta;
      }

      void diffuse ( double dt )
      {
        for ( int i = 0; i < m; ++i )
          for ( int j = 0; j < n; ++j )
            {
              Node  *pn = getNode(i,j);
              if ( pn->b&2 )
                {
                  pn->h += 0.5*dt*(pn->u+pn->qp);
                  pn->l = pn->h;
                }
            }
        for ( int i = 0; i < m; ++i )
          for ( int j = 0; j < n; ++j )
            {
              Node  *pn = getNode(i,j);
              if ( pn->b&2 )
                {
                  if ( getNodeP(i+1,j)->b&2 )
                    {
                      double  dh = dt*odiff*(pn->l-getNodeP(i+1,j)->l);
                      pn->h -= dh;
                      if ( getNodeP(i+1,j)->b==2 ) getNodeP(i+1,j)->h += dh;
                    }
                  if ( getNodeP(i,j+1)->b&2 )
                    {
                      double  dh = dt*odiff*(pn->l-getNodeP(i,j+1)->l);
                      pn->h -= dh;
                      if ( getNodeP(i,j+1)->b==2 )  getNodeP(i,j+1)->h += dh;
                    }
                  pn->h += 0.5*dt*(pn->u+pn->qp);
                }
            }
      }
  };

  class Connector
  {
    public:
      double  hscale, x0, y0, alpha, dx0, dy0, dalpha;
      vector<vector<double>>  x, y, vx, vy, vz;
      OceanGrid  *g;

      Connector ( OceanGrid *g, double hscale = 1, double x0 = 0, double y0 = 0, double alpha = 0 )
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

      void convertVelocities(bool minimize_advection = true)
      {

        //todo: add and substract 4 of b
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

              AssertThrow(!std::isnan(vx[i][j]) && !std::isnan(vy[i][j]),
                          aspect::ExcMessage("vx and/or vy are nan at " + std::to_string(i) + ":" + std::to_string(j) + ", vxx:y = "  + std::to_string(vx[i][j]) + ":" + std::to_string(vy[i][j])));
// Only continental area
              if ( g->getNode(i,j)->b == 0 )
                {
                  ++n;
                  xmean += i;
                  ymean += j;
                  vxmean += vx[i][j];
                  vymean += vy[i][j];
                  cross += i*vy[i][j]-j*vx[i][j];
                  //if (std::fabs(vx[i][j]-(-y[i][j])) > std::numeric_limits<double>::epsilon() || std::fabs(vy[i][j]-x[i][j]) > std::numeric_limits<double>::epsilon())
                  //  {
                  //    std::cout << "wrong: i:j = " + std::to_string(i) + ":" + std::to_string(j) + ", x:y = " << x[i][j] << ":" << y[i][j] << ", vx:vy = " <<vx[i][j] << ":" << vy[i][j] << ", diff x:y = "  << vx[i][j]-(-y[i][j]) << ":" << vy[i][j]-x[i][j] << std::endl;
                  //    AssertThrow(false, aspect::ExcMessage("wrong!!"));
                  //  }
                  rsq += i*i+j*j;
                }
            }
        if (n == 0 || minimize_advection == false)
          {
            return;
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
        //for ( int i = 0; i < g->m; i=i+10 )
        //  for ( int j = 0; j < g->n; j=j+10 )
        //    {
        //      if (
        //        //(i == 5 && j == 5) ||
        //        //(i == 50 && j == 50) ||
        //        //(i == 100 && j == 100) ||
        //        (i == 116 && j == 26)
        //      )
        //        {
        //          std::cout << " vx:vy = " << i << ":" << j << " = "<< vx[i][j] << " : " << vy[i][j] << std::endl;
        //        }
        //    }
        for ( int i = 0; i < g->m; ++i )
          for ( int j = 0; j < g->n; ++j )
            {
              //if (
              //  //(i == 5 && j == 5) ||
              //  //(i == 50 && j == 50) ||
              //  //(i == 100 && j == 100) ||
              //  (i == 116 && j == 26)
              //)
              //  {
              //    std::cout << " i:j= " << i << ":" << j << " vx:vy = "<< vx[i][j] << " : " << vy[i][j] << ", vx-dx0-dalpha*j = " << vx[i][j]-(dx0-dalpha *j) << ", vy-dy0+daplha*i = " << vy[i][j]-(dy0+dalpha *i) << std::endl;
              //  }
              vx[i][j] -= dx0-dalpha*j;
              vy[i][j] -= dy0+dalpha*i;
              AssertThrow(!std::isnan(vx[i][j]) && !std::isnan(vy[i][j]),
                          aspect::ExcMessage("vx and/or vy are nan at " + std::to_string(i) + ":" + std::to_string(j) + ", vxx:y = "  + std::to_string(vx[i][j]) + ":" + std::to_string(vy[i][j])));
            }

        double dx0_old = dx0;
        dx0 = (c*dx0-s*dy0)*hscale;
        dy0 = (s*dx0_old+c*dy0)*hscale;

        //std::cout << "residual x:y = " << 0 << ":" << 0 << " = "<< vx[0][0] i                     << " : " << vy[0][0] << std::endl;
        //std::cout << "residual x:y = " << 0 << ":" << openlem_ny << " = "<< vx[0][openlem_ny-1]            << " : " << vy[0][openlem_ny-1] << std::endl;
        //std::cout << "residual x:y = " << 0 << ":" << 0 << " = "<< vx[openlem_nx-1][openlem_ny-1] << " : " << vy[openlem_nx-1][openlem_ny-1] << std::endl;
        //std::cout << "residual x:y = " << 0 << ":" << 0 << " = "<< vx[openlem_nx-1][0]            << " : " << vy[openlem_nx-1][0] << std::endl;
        //for ( int i = 0; i < g->m; i=i+10 )
        //  for ( int j = 0; j < g->n; j=j+10 )
        //    {
        //      //if (vx[i][j] > 0 && vy[i][j] > 0)
        //      if (
        //        //(i == 5 && j == 5) ||
        //        //(i == 50 && j == 50) ||
        //        //(i == 100 && j == 100) ||
        //        (i == 116 && j == 126)
        //      )
        //        {
        //          std::cout << "residual x:y = " << i << ":" << j << " = "<< vx[i][j] << " : " << vy[i][j] << std::endl;
        //        }
        //    }
      }

      void computeUpliftRate()
      {
        std::cout << "connector computeUpliftRate" << std::endl;
        for ( int i = 0; i < g->m; ++i )
          for ( int j = 0; j < g->n; ++j )
            if ( g->getNode(i,j)->b != 3 ) //&1 )
              //if ( !g->getNode(i,j)->b &1 )
              {
                double  dhdx = vx[i][j]<0 ?
                               g->getNodeP(i+1,j)->h-g->getNode(i,j)->h :
                               g->getNode(i,j)->h-g->getNodeP(i-1,j)->h;
                double  dhdy = vy[i][j]<0 ?
                               g->getNodeP(i,j+1)->h-g->getNode(i,j)->h :
                               g->getNode(i,j)->h-g->getNodeP(i,j-1)->h;
                //if (i == 19 && j == 24)
                //  std::cout << "before dhdx:y = " << dhdx << ':' << dhdy << ", u = " << g->getNode(i,j)->u << ", vx:y:z = " << vx[i][j] << ":" << vy[i][j] << ":" << vz[i][j] << std::endl;
                g->getNode(i,j)->u = vz[i][j]-vx[i][j]*dhdx-vy[i][j]*dhdy;
                //if (i == 19 && j == 24)
                //  std::cout << "dhdx:y = " << dhdx << ':' << dhdy << ", u = " << g->getNode(i,j)->u << ", vx:y:z = " << vx[i][j] << ":" << vy[i][j] << ":" << vz[i][j] << std::endl;
              }
      }

      void updateCoordinateSystem ( double dt )
      {
        std::cout << ".. connector updateCoordinatesystem ..";
        //std::cout << "connector updateCoordinatesystem x0:y0 = " << x0 << ":" << y0 << ", dx0:dy0 = " << dx0 << ":" << dy0 << ", dt = " << dt << ", aplha = " << alpha << ", dalhpa = " << dalpha <<std::endl;
        x0 += dx0*dt;
        y0 += dy0*dt;
        alpha += dalpha*dt;
        //std::cout << "connector updateCoordinatesystem x0:y0 = " << x0 << ":" << y0 << ", dx0:dy0 = " << dx0 << ":" << dy0 << ", dt = " << dt << ", alpha = " << alpha << std::endl;
        updateXY();
      }

      double interpolation(double x, double y, const std::vector<std::vector<double>> &interperolat, const std::vector<std::vector<double>> &interperolat2,  double dx, double dy, double origin_x, double origin_y) const
      {
        /*
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
          }*/
        // find the closest point
        double closest_distance = std::numeric_limits<double>::infinity();
        size_t closest_xi = std::numeric_limits<size_t>::signaling_NaN();
        size_t closest_yi = std::numeric_limits<size_t>::signaling_NaN();
        for (unsigned int xi = 0; xi < interperolat.size()-1; ++xi)
          {
            for (unsigned int yi = 0; yi < interperolat[xi].size()-1; ++yi)
              {
                double distance = (x-xi*dx-origin_x)*(x-xi*dx-origin_x)+(y-yi*dy-origin_y)*(y-yi*dy-origin_y);
                if (distance < closest_distance)
                  {
                    closest_distance = distance;
                    closest_xi = xi;
                    closest_yi = yi;
                  }

              }
          }

        //if (x > -6345.87-100. && x < -6345.87+4000. && y > 28103. - 200. && y < 28103. + 200. )
        //if (x > -6345.87-100. && x < -6345.87+5000. && y > 28103. - 400. && y < 28103. + 400. )
        //if (x > 3125-100. && x < 3125+5000. && y > 45312.5 - 400. && y < 45312.5 + 400. )
        //std::cout << "closest_xi = " << closest_xi << ", closest_yi = " << closest_yi << ", interp = " << interperolat[closest_xi][closest_yi] << ", interp2 = " <<  interperolat2[closest_xi][closest_yi] << std::endl;
        return interperolat[closest_xi][closest_yi];// - interperolat2[closest_xi][closest_yi];//g->getNode(i,j)->h;//interperolat[i][j];
      }

      void write_vtk(double reference_surface_height, int timestep, double time, std::string path, std::string prestring = "") const
      {
        std::string timestep_string = std::to_string(timestep);
        size_t n_zero = 5;
        std::string relative_path = "openLEM/openlem_surface" + prestring + std::string(n_zero - std::min(n_zero, timestep_string.length()), '0') +timestep_string +".vtu";
        std::string full_path = path + relative_path;
        std::string pvtu_full_path = path + "openlem_surface" + prestring + ".pvd";
        std::cout << "writing openlem VTK file to " << full_path << std::endl;
        std::vector<double> grid_x(0);
        std::vector<double> grid_y(0);
        std::vector<double> grid_z(0);
        std::vector<double> grid_elevation(0);
        std::vector<double> grid_uplift(0);
        std::vector<double> grid_boundary(0);
        std::vector<double> grid_q(0);
        std::vector<double> grid_vx(0);
        std::vector<double> grid_vy(0);

        std::vector<std::vector<size_t>> grid_connectivity(0);

        unsigned int n_cell = x.size() * y.size();
        unsigned int n_p = (x.size() + 1) * (y.size() + 1);


        grid_x.resize(n_p);
        grid_z.resize(n_p);
        grid_y.resize(n_p);

        grid_elevation.resize(n_p);
        grid_uplift.resize(n_p);
        grid_boundary.resize(n_p);
        grid_q.resize(n_p);
        grid_vx.resize(n_p);
        grid_vy.resize(n_p);
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
                grid_uplift[counter] = g->getNode(i,j)->u;
                grid_boundary[counter] = g->getNode(i,j)->b;
                grid_q[counter] = g->getNode(i,j)->q;
                grid_vx[counter] = vx[i][j];//g->getNode(i,j)->u;
                grid_vy[counter] = vy[i][j];//g->getNode(i,j)->u;
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
        myfile.open (full_path);
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
        buffer << R"(<DataArray type="Float32" Name="Uplift" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer <<  grid_uplift[i] << std::endl;
          }
        buffer << "</DataArray>" << std::endl;
        buffer << R"(<DataArray type="Float32" Name="Boundary value" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer <<  grid_boundary[i] << std::endl;
          }
        buffer << "</DataArray>" << std::endl;
        buffer << R"(<DataArray type="Float32" Name="q" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer <<  grid_q[i] << std::endl;
          }
        buffer << "</DataArray>" << std::endl;
        buffer << R"(<DataArray type="Float32" Name="vx" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer <<  grid_vx[i] << std::endl;
          }
        buffer << "</DataArray>" << std::endl;
        buffer << R"(<DataArray type="Float32" Name="vy" format="ascii">)" << std::endl;
        for (size_t i = 0; i < n_p; ++i)
          {
            buffer <<  grid_vy[i] << std::endl;
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

        myfile.close();
        if (prestring != "")
          {
            std::cout << "Finished writing openlem VTK file" << std::endl;
            return;
          }
        std::string new_pvtu_line = "    <DataSet timestep=\"" + std::to_string(time) + "\" group=\"\" part=\"0\" file=\"" + relative_path + "\"/>";
        if (timestep == 0)
          {

            myfile.open (pvtu_full_path);

            buffer << R"(<?xml version="1.0"?>)" << std::endl;
            buffer << R"(<VTKFile type="Collection" version="0.1" ByteOrder="LittleEndian">)";
            buffer << R"(  <Collection>)" << std::endl;
            buffer << new_pvtu_line << std::endl;
            buffer << R"(  </Collection>)" << std::endl;
            buffer << R"(</VTKFile>)" << std::endl;
            myfile << buffer.str();
            buffer.str(std::string());
            myfile.close();
          }
        else
          {
            // read file into vector
            std::ifstream infile(pvtu_full_path);
            std::string line;
            std::vector<std::string> lines_vector;
            while (std::getline(infile, line))
              {
                std::istringstream iss(line);
                lines_vector.emplace_back(line);
              }
            infile.close();
            lines_vector.insert(lines_vector.end() - 2, new_pvtu_line + '\n');

            myfile.open (pvtu_full_path);

            for (std::string &line : lines_vector)
              {
                buffer << line << std::endl;
              }
            myfile << buffer.str();
            buffer.str(std::string());
          }

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
        void execute_openlem(openlem::OceanGrid &grid,
                             //std::vector<double> &elevation,
                             //std::vector<double> &extra_vtk_field,
                             //std::vector<double> &velocity_x,
                             //std::vector<double> &velocity_y,
                             //std::vector<double> &velocity_z,
                             const double &openlem_timestep_in_years,
                             const unsigned int &openlem_iterations,
                             const std::string &dirname);

        /**
         * Function to fill the openlem arrays (height and velocities) with the data received from ASPECT in the correct index order.
         */
        void fill_openlem_arrays(openlem::OceanGrid &grid,
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
        openlem::OceanGrid grid_old;
        openlem::OceanGrid grid_new;
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
        openlem::Point deepest_point;
        double openlem_ocean_diffusivity;

        bool openlem_minimize_advection;
        double openlem_kd;
        double openlem_kt;
        //std::vector<std::vector<double>> mesh_velocity_z;
        std::vector<std::vector<double>> aspect_mesh_dh;
        std::vector<std::vector<double>> aspect_mesh_z;
        std::vector<std::vector<typename DoFHandler<dim>::active_cell_iterator>> aspect_mesh_cell_id;
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
