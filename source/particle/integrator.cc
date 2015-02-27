/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n).
       * This requires only one step per integration, and doesn't involve any extra data.
       */
      template <int dim, class T>
      class EulerIntegrator : public Interface<dim, T>
      {
        public:
          virtual bool integrate_step(Particle::World<dim, T> *world, const double dt)
          {
            typename std::multimap<LevelInd, T> &particles = world->get_particles();
            typename std::multimap<LevelInd, T>::iterator       it;
            Point<dim>                          loc, vel;

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                it->second.set_location(loc + dt*vel);
              }

            return false;
          };
          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info) {};
          virtual unsigned int data_len() const
          {
            return 0;
          };
          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            return pos;
          };
          virtual void write_data(std::vector<double> &data, const double &id_num) const
          {
          };
      };

      /**
       * Runge Kutta second order integrator, where y_{n+1} = y_n + dt*v(0.5*k_1), k_1 = dt*v(y_n).
       * This scheme requires storing the original location, and the read/write_data functions reflect this.
       */
      template <int dim, class T>
      class RK2Integrator : public Interface<dim, T>
      {
        private:
          unsigned int                    step;
          std::map<double, Point<dim> >   loc0;

        public:
          RK2Integrator(void)
          {
            step = 0;
            loc0.clear();
          };

          virtual bool integrate_step(Particle::World<dim, T> *world, const double dt)
          {
            typename std::multimap<LevelInd, T> &particles = world->get_particles();
            typename std::multimap<LevelInd, T>::iterator       it;
            Point<dim>                          loc, vel;
            double                              id_num;

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                id_num = it->second.get_id();
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                if (step == 0)
                  {
                    loc0[id_num] = loc;
                    it->second.set_location(loc + 0.5*dt*vel);
                  }
                else if (step == 1)
                  {
                    it->second.set_location(loc0[id_num] + dt*vel);
                  }
                else
                  {
                    // Error!
                  }
              }

            if (step == 1) loc0.clear();
            step = (step+1)%2;

            // Continue until we're at the last step
            return (step != 0);
          };

          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {
            // Add the loc0 data
            data_info.push_back(MPIDataInfo("loc0", dim));
          };

          virtual unsigned int data_len() const
          {
            return dim;
          };

          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            unsigned int    i, p = pos;

            // Read location data
            for (i=0; i<dim; ++i)
              {
                loc0[id_num](i) = data[p++];
              }

            return p;
          };

          virtual void write_data(std::vector<double> &data, const double &id_num) const
          {
            unsigned int    i;
            typename std::map<double, Point<dim> >::const_iterator it;

            // Write location data
            it = loc0.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
          };
      };


      /**
       * Runge Kutta fourth order integrator, where y_{n+1} = y_n + (1/6)*k1 + (1/3)*k2 + (1/3)*k3 + (1/6)*k4
       * and k1, k2, k3, k4 are defined as usual.
       * This scheme requires storing the original location and intermediate k1, k2, k3 values,
       * so the read/write_data functions reflect this.
       */
      template <int dim, class T>
      class RK4Integrator : public Interface<dim, T>
      {
        private:
          unsigned int                    step;
          std::map<double, Point<dim> >   loc0, k1, k2, k3;

        public:
          RK4Integrator(void)
          {
            step = 0;
            loc0.clear();
            k1.clear();
            k2.clear();
            k3.clear();
          };

          virtual bool integrate_step(Particle::World<dim, T> *world, const double dt)
          {
            typename std::multimap<LevelInd, T> &particles = world->get_particles();
            typename std::multimap<LevelInd, T>::iterator       it;
            Point<dim>                          loc, vel, k4;
            double                              id_num;

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                id_num = it->second.get_id();
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                if (step == 0)
                  {
                    loc0[id_num] = loc;
                    k1[id_num] = dt*vel;
                    it->second.set_location(loc + 0.5*k1[id_num]);
                  }
                else if (step == 1)
                  {
                    k2[id_num] = dt*vel;
                    it->second.set_location(loc0[id_num] + 0.5*k2[id_num]);
                  }
                else if (step == 2)
                  {
                    k3[id_num] = dt*vel;
                    it->second.set_location(loc0[id_num] + k3[id_num]);
                  }
                else if (step == 3)
                  {
                    k4 = dt*vel;
                    it->second.set_location(loc0[id_num] + (k1[id_num]+2.0*k2[id_num]+2.0*k3[id_num]+k4)/6.0);
                  }
                else
                  {
                    // Error!
                  }
              }

            step = (step+1)%4;
            if (step == 0)
              {
                loc0.clear();
                k1.clear();
                k2.clear();
                k3.clear();
              }

            // Continue until we're at the last step
            return (step != 0);
          };

          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {
            // Add the loc0, k1, k2, and k3 data
            data_info.push_back(MPIDataInfo("loc0", dim));
            data_info.push_back(MPIDataInfo("k1", dim));
            data_info.push_back(MPIDataInfo("k2", dim));
            data_info.push_back(MPIDataInfo("k3", dim));
          };

          virtual unsigned int data_len() const
          {
            return 4*dim;
          };

          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            unsigned int    i, p = pos;

            // Read location data
            for (i=0; i<dim; ++i)
              {
                loc0[id_num](i) = data[p++];
              }
            // Read k1, k2 and k3
            for (i=0; i<dim; ++i)
              {
                k1[id_num](i) = data[p++];
              }
            for (i=0; i<dim; ++i)
              {
                k2[id_num](i) = data[p++];
              }
            for (i=0; i<dim; ++i)
              {
                k3[id_num](i) = data[p++];
              }

            return p;
          };

          virtual void write_data(std::vector<double> &data, const double &id_num) const
          {
            typename std::map<double, Point<dim> >::const_iterator it;
            unsigned int  i;

            // Write location data
            it = loc0.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            // Write k1, k2 and k3
            it = k1.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            it = k2.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            it = k3.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
          };
      };

      /**
       * Integrator which chooses Euler, RK2 or RK4 depending on characteristics of the cell a particle is in.
       * Currently used for research only.
       */
      template <int dim, class T>
      class HybridIntegrator : public Interface<dim, T>
      {
        private:
          enum IntegrationScheme
          {
            SCHEME_UNDEFINED,
            SCHEME_EULER,
            SCHEME_RK2,
            SCHEME_RK4
          };

          unsigned int                    step;
          std::map<double, Point<dim> >   loc0, k1, k2, k3;
          std::map<double, IntegrationScheme>        scheme;

          virtual IntegrationScheme select_scheme(const std::vector<Point<dim> > &cell_vertices, const std::vector<Point<dim> > &cell_velocities, const double timestep)
          {
            return cell_vertices[0][0] > 0.5 ? SCHEME_RK4 : SCHEME_EULER;
          };

        public:
          HybridIntegrator()
          {
            step = 0;
            loc0.clear();
            k1.clear();
            k2.clear();
            k3.clear();
            scheme.clear();
          };

          virtual bool integrate_step(Particle::World<dim, T> *world, const double dt)
          {
            typename std::multimap<LevelInd, T> &particles = world->get_particles();
            typename std::multimap<LevelInd, T>::iterator    it;
            const DoFHandler<dim>                            *dh = world->get_dof_handler();
            const Mapping<dim>                               *mapping = world->get_mapping();
            const parallel::distributed::Triangulation<dim>  *tria = world->get_triangulation();
            const LinearAlgebra::BlockVector         *solution = world->get_solution();
            Point<dim>                                       loc, vel, k4;
            double                                           id_num;
            LevelInd                                         cur_level_ind;
            IntegrationScheme                                cur_scheme;
            typename DoFHandler<dim>::active_cell_iterator   found_cell;
            Functions::FEFieldFunction<dim, DoFHandler<dim>, LinearAlgebra::BlockVector> fe_value(*dh, *solution, *mapping);

            // If this is the first step, go through all the cells and determine
            // which integration scheme the particles in each cell should use
            if (step == 0)
              {
                Vector<double>                single_res(dim+2);
                std::vector<Point<dim> >    cell_vertices, cell_velocities;
                std::vector<Vector<double> > temp_vals;

                cell_vertices.resize(GeometryInfo<dim>::vertices_per_cell);
                cell_velocities.resize(GeometryInfo<dim>::vertices_per_cell);
                cur_level_ind.first = cur_level_ind.second = -1;
                cur_scheme = SCHEME_UNDEFINED;
                temp_vals.resize(GeometryInfo<dim>::vertices_per_cell, single_res);

                for (it=particles.begin(); it!=particles.end(); ++it)
                  {
                    if (cur_level_ind != it->first)
                      {
                        cur_level_ind = it->first;
                        found_cell = typename DoFHandler<dim>::active_cell_iterator(tria, cur_level_ind.first, cur_level_ind.second, dh);

                        // Ideally we should use all quadrature point velocities, but for now
                        // we just evaluate them at the vertices
                        fe_value.set_active_cell(found_cell);
                        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i) cell_vertices[i] = found_cell->vertex(i);
                        fe_value.vector_value_list(cell_vertices, temp_vals);
                        for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
                          {
                            for (unsigned int d=0; d<dim; ++d)
                              {
                                cell_velocities[i][d] = temp_vals[i][d];
                              }
                          }

                        cur_scheme = select_scheme(cell_vertices, cell_velocities, dt);
                      }
                    scheme[it->second.get_id()] = cur_scheme;
                  }
              }

            for (it=particles.begin(); it!=particles.end(); ++it)
              {
                id_num = it->second.get_id();
                loc = it->second.get_location();
                vel = it->second.get_velocity();
                switch (scheme[id_num])
                  {
                    case SCHEME_EULER:
                      if (step == 0)
                        {
                          it->second.set_location(loc + dt*vel);
                        }
                      it->second.set_vel_check(false);
                      break;
                    case SCHEME_RK2:
                      if (step == 0)
                        {
                          loc0[id_num] = loc;
                          it->second.set_location(loc + 0.5*dt*vel);
                        }
                      else if (step == 1)
                        {
                          it->second.set_location(loc0[id_num] + dt*vel);
                        }
                      if (step != 0) it->second.set_vel_check(false);
                      break;
                    case SCHEME_RK4:
                      if (step == 0)
                        {
                          loc0[id_num] = loc;
                          k1[id_num] = dt*vel;
                          it->second.set_location(loc + 0.5*k1[id_num]);
                        }
                      else if (step == 1)
                        {
                          k2[id_num] = dt*vel;
                          it->second.set_location(loc0[id_num] + 0.5*k2[id_num]);
                        }
                      else if (step == 2)
                        {
                          k3[id_num] = dt*vel;
                          it->second.set_location(loc0[id_num] + k3[id_num]);
                        }
                      else if (step == 3)
                        {
                          k4 = dt*vel;
                          it->second.set_location(loc0[id_num] + (k1[id_num]+2.0*k2[id_num]+2.0*k3[id_num]+k4)/6.0);
                        }
                      break;
                    default:
                      AssertThrow(false, ExcMessage("Unknown integration scheme for hybrid integrator."));
                      break;
                  }
              }

            step = (step+1)%4;
            if (step == 0)
              {
                loc0.clear();
                k1.clear();
                k2.clear();
                k3.clear();
                scheme.clear();
              }

            // Continue until we're at the last step
            return (step != 0);
          };

          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {
            // Add the loc0, k1, k2, and k3 data
            data_info.push_back(MPIDataInfo("loc0", dim));
            data_info.push_back(MPIDataInfo("k1", dim));
            data_info.push_back(MPIDataInfo("k2", dim));
            data_info.push_back(MPIDataInfo("k3", dim));
            data_info.push_back(MPIDataInfo("scheme", 1));
          };

          virtual unsigned int data_len() const
          {
            return (4*dim+1);
          };

          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
          {
            unsigned int    i, p = pos;

            // Read location data
            for (i=0; i<dim; ++i)
              {
                loc0[id_num](i) = data[p++];
              }
            // Read k1, k2 and k3
            for (i=0; i<dim; ++i)
              {
                k1[id_num](i) = data[p++];
              }
            for (i=0; i<dim; ++i)
              {
                k2[id_num](i) = data[p++];
              }
            for (i=0; i<dim; ++i)
              {
                k3[id_num](i) = data[p++];
              }
            scheme[id_num] = (IntegrationScheme)data[p++];

            return p;
          };

          virtual void write_data(std::vector<double> &data, const double &id_num) const
          {
            typename std::map<double, Point<dim> >::const_iterator it;
            typename std::map<double, IntegrationScheme>::const_iterator sit;
            unsigned int  i;

            // Write location data
            it = loc0.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            // Write k1, k2 and k3
            it = k1.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            it = k2.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            it = k3.find(id_num);
            for (i=0; i<dim; ++i)
              {
                data.push_back(it->second(i));
              }
            // Write the integration scheme for this particle
            sit = scheme.find(id_num);
            data.push_back((double)(sit->second));
          };
      };

      template <int dim, class T>
      Interface<dim,T> *
      create_integrator_object (const std::string &integrator_name)
      {
        if (integrator_name == "euler")
          return new EulerIntegrator<dim,T>();
        else if (integrator_name == "rk2")
          return new RK2Integrator<dim,T>();
        else if (integrator_name == "rk4")
          return new RK4Integrator<dim,T>();
        else if (integrator_name == "hybrid")
          return new HybridIntegrator<dim,T>();
        else
          Assert (false, ExcNotImplemented());

        return 0;
      }


      std::string
      integrator_object_names ()
      {
        return ("euler|"
                "rk2|"
                "rk4|"
                "hybrid");
      }


      // explicit instantiations
      template
      Interface<2,Particle::BaseParticle<2> > *
      create_integrator_object (const std::string &integrator_name);
      template
      Interface<3,Particle::BaseParticle<3> > *
      create_integrator_object (const std::string &integrator_name);
    }
  }
}
