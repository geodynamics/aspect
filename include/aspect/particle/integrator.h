/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id$  */

#ifndef __aspect__particle_integrator_h
#define __aspect__particle_integrator_h

#include <aspect/particle/particle.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    /**
     * Integrator is an abstract class defining virtual methods for performing integration of
     * particle paths through the simulation velocity field.
     */
    template <int dim, class T>
    class Integrator
    {
      public:
        virtual ~Integrator(void) {};

        /**
         * Perform an integration step of moving the particles by the specified timestep dt.
         * Implementations of this function must update the particle location.
         * If the integrator requires multiple internal steps, this function must return
         * true until all internal steps are finished. Between calls to this function
         * the velocity at the updated particle positions is evaluated and passed to
         * integrate_step during the next call.
         *
         * @param [in,out] particles The set of particles to integrate. The particle positions
         * will be changed in this function based on the integration scheme.
         * @param [in] dt The timestep length to perform the integration.
         * @return Whether this function needs to be called again (true) for additional
         * integration steps or if all internal steps are complete (false).
         */
        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, const double dt) = 0;

        /**
         * Secify the MPI types and data sizes involved in transferring integration related
         * information between processes. If the integrator samples velocities at different
         * locations and the particle moves between processes during the integration step,
         * the sampled velocities must be transferred with the particle.
         *
         * @param [in,out] data_info Adds MPI data info to the specified vector indicating
         * the quantity and type of values the integrator needs saved for this particle.
         */
        virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info) = 0;

        /**
        * Return data length of the integration related data required for communication
        * in terms of number of doubles.
        *
        * @return The number of doubles required to store the relevant integrator data.
        */
        virtual unsigned int data_len() const = 0;

        /**
         * Read integration related data for a particle specified by id_num from the data vector.
         *
         * @param [in] data The vector of double data to read from.
         * @param [in] pos The position in the data vector to start reading from.
         * @param [in] id_num The id number of the particle to read the data for.
         * @return The position in the vector of the next unread double.
         */
        virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num) = 0;

        /**
         * Write integration related data to a vector for a particle
         * specified by id_num.
         *
         * @param [in,out] data The vector of doubles to write integrator data into.
         * @param [in] id_num The id number of the particle to read the data for.
         */
        virtual void write_data(std::vector<double> &data, const double &id_num) const = 0;
    };

    /**
     * Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n).
     * This requires only one step per integration, and doesn't involve any extra data.
     */
    template <int dim, class T>
    class EulerIntegrator : public Integrator<dim, T>
    {
      public:
        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, const double dt)
        {
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
    class RK2Integrator : public Integrator<dim, T>
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

        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, const double dt)
        {
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
    class RK4Integrator : public Integrator<dim, T>
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

        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, const double dt)
        {
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
                  it->second.set_location(loc0[id_num] + (k1[id_num]+2*k2[id_num]+2*k3[id_num]+k4)/6.0);
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
    class HybridIntegrator : public Integrator<dim, T>
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
        const parallel::distributed::Triangulation<dim>   *tria;
        const DoFHandler<dim>           *dh;
        const Mapping<dim>              *mapping;
        const TrilinosWrappers::MPI::BlockVector *solution;

        virtual IntegrationScheme select_scheme(const std::vector<Point<dim> > &cell_vertices, const std::vector<Point<dim> > &cell_velocities, const double timestep)
        {
          return cell_vertices[0][0] > 0.5 ? SCHEME_RK4 : SCHEME_EULER;
        };

      public:
        HybridIntegrator(const parallel::distributed::Triangulation<dim> *new_tria, const DoFHandler<dim> *new_dh, const Mapping<dim> *new_mapping, const TrilinosWrappers::MPI::BlockVector *new_solution)
        {
          step = 0;
          loc0.clear();
          k1.clear();
          k2.clear();
          k3.clear();
          scheme.clear();
          tria = new_tria;
          dh = new_dh;
          mapping = new_mapping;
          solution = new_solution;
        };

        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, const double dt)
        {
          typename std::multimap<LevelInd, T>::iterator       it;
          Point<dim>                          loc, vel, k4;
          double                              id_num;
          LevelInd                            cur_level_ind;
          IntegrationScheme                   cur_scheme;
          typename DoFHandler<dim>::active_cell_iterator found_cell;
          Functions::FEFieldFunction<dim, DoFHandler<dim>, TrilinosWrappers::MPI::BlockVector> fe_value(*dh, *solution, *mapping);

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
                        it->second.set_location(loc0[id_num] + (k1[id_num]+2*k2[id_num]+2*k3[id_num]+k4)/6.0);
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
  }
}

#endif

