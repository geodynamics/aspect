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
    // Integrator is an abstract class defining virtual methods for performing integration of
    // particle paths through the simulation velocity field
    template <int dim, class T>
    class Integrator
    {
      public:
        virtual ~Integrator(void) {};

        // Perform an integration step of moving the particles by the specified timestep dt.
        // Implementations of this function should update the particle location.
        // If the integrator requires multiple internal steps, this should return true until
        // all internal steps are finished.  Between calls to this, the particles will be moved
        // based on the velocities at their newly specified positions.
        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, double dt) = 0;

        // Specify the MPI types and data sizes involved in transferring integration related
        // information between processes. If the integrator samples velocities at different
        // locations and the particle moves between processes during the integration step,
        // the sampled velocities must be transferred with the particle.
        virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info) = 0;

        // Must return data length in bytes of the integration related data
        // for a particle in a given format.
        virtual unsigned int data_len(ParticleDataFormat format) const = 0;

        // Read integration related data for a particle specified by id_num
        // Returns the data pointer updated to point to the next unwritten byte
        virtual const char *read_data(ParticleDataFormat format, double id_num, const char *data) = 0;

        // Write integration related data for a particle specified by id_num
        // Returns the data pointer updated to point to the next unwritten byte
        virtual char *write_data(ParticleDataFormat format, double id_num, char *data) = 0;
    };

    // Euler scheme integrator, where y_{n+1} = y_n + dt * v(y_n)
    // This requires only one step per integration, and doesn't involve any extra data.
    template <int dim, class T>
    class EulerIntegrator : public Integrator<dim, T>
    {
      public:
        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, double dt)
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
        virtual unsigned int data_len(ParticleDataFormat format) const
        {
          return 0;
        };
        virtual const char *read_data(ParticleDataFormat format, double id_num, const char *data)
        {
          return data;
        };
        virtual char *write_data(ParticleDataFormat format, double id_num, char *data)
        {
          return data;
        };
    };

    // Runge Kutta second order integrator, where y_{n+1} = y_n + dt*f(k_1, k_1 = v(
    // This scheme requires storing the original location, and the read/write_data functions reflect this
    template <int dim, class T>
    class RK2Integrator : public Integrator<dim, T>
    {
      private:
        unsigned int                    _step;
        std::map<double, Point<dim> >   _loc0;

      public:
        RK2Integrator(void)
        {
          _step = 0;
          _loc0.clear();
        };

        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, double dt)
        {
          typename std::multimap<LevelInd, T>::iterator       it;
          Point<dim>                          loc, vel;
          double                              id_num;

          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              id_num = it->second.get_id();
              loc = it->second.get_location();
              vel = it->second.get_velocity();
              if (_step == 0)
                {
                  _loc0[id_num] = loc;
                  it->second.set_location(loc + 0.5*dt*vel);
                }
              else if (_step == 1)
                {
                  it->second.set_location(_loc0[id_num] + dt*vel);
                }
              else
                {
                  // Error!
                }
            }

          if (_step == 1) _loc0.clear();
          _step = (_step+1)%2;

          // Continue until we're at the last step
          return (_step != 0);
        };

        virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Add the _loc0 data
          data_info.push_back(MPIDataInfo("loc0", dim, MPI_DOUBLE, sizeof(double)));
        };

        virtual unsigned int data_len(ParticleDataFormat format) const
        {
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                return dim*sizeof(double);
            }
          return 0;
        };

        virtual const char *read_data(ParticleDataFormat format, double id_num, const char *data)
        {
          const char     *p = data;
          unsigned int    i;

          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Read location data
                for (i=0; i<dim; ++i)
                  {
                    _loc0[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                break;
            }

          return p;
        };

        virtual char *write_data(ParticleDataFormat format, double id_num, char *data)
        {
          char          *p = data;
          unsigned int  i;

          // Then write our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Write location data
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _loc0[id_num](i);
                    p += sizeof(double);
                  }
                break;
            }

          return p;
        };
    };


    template <int dim, class T>
    class RK4Integrator : public Integrator<dim, T>
    {
      private:
        unsigned int                    _step;
        std::map<double, Point<dim> >   _loc0, _k1, _k2, _k3;

      public:
        RK4Integrator(void)
        {
          _step = 0;
          _loc0.clear();
          _k1.clear();
          _k2.clear();
          _k3.clear();
        };

        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, double dt)
        {
          typename std::multimap<LevelInd, T>::iterator       it;
          Point<dim>                          loc, vel, k4;
          double                              id_num;

          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              id_num = it->second.get_id();
              loc = it->second.get_location();
              vel = it->second.get_velocity();
              if (_step == 0)
                {
                  _loc0[id_num] = loc;
                  _k1[id_num] = dt*vel;
                  it->second.set_location(loc + 0.5*_k1[id_num]);
                }
              else if (_step == 1)
                {
                  _k2[id_num] = dt*vel;
                  it->second.set_location(_loc0[id_num] + 0.5*_k2[id_num]);
                }
              else if (_step == 2)
                {
                  _k3[id_num] = dt*vel;
                  it->second.set_location(_loc0[id_num] + _k3[id_num]);
                }
              else if (_step == 3)
                {
                  k4 = dt*vel;
                  it->second.set_location(_loc0[id_num] + (_k1[id_num]+2*_k2[id_num]+2*_k3[id_num]+k4)/6.0);
                }
              else
                {
                  // Error!
                }
            }

          _step = (_step+1)%4;
          if (_step == 0)
            {
              _loc0.clear();
              _k1.clear();
              _k2.clear();
              _k3.clear();
            }

          // Continue until we're at the last step
          return (_step != 0);
        };

        virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Add the _loc0, _k1, _k2, and _k3 data
          data_info.push_back(MPIDataInfo("loc0", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("k1", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("k2", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("k3", dim, MPI_DOUBLE, sizeof(double)));
        };

        virtual unsigned int data_len(ParticleDataFormat format) const
        {
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                return 4*dim*sizeof(double);
            }
          return 0;
        };

        virtual const char *read_data(ParticleDataFormat format, double id_num, const char *data)
        {
          const char     *p = data;
          unsigned int    i;

          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Read location data
                for (i=0; i<dim; ++i)
                  {
                    _loc0[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                // Read k1, k2 and k3
                for (i=0; i<dim; ++i)
                  {
                    _k1[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    _k2[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    _k3[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                break;
            }

          return p;
        };

        virtual char *write_data(ParticleDataFormat format, double id_num, char *data)
        {
          char          *p = data;
          unsigned int  i;

          // Then write our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Write location data
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _loc0[id_num](i);
                    p += sizeof(double);
                  }
                // Write k1, k2 and k3
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _k1[id_num](i);
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _k2[id_num](i);
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _k3[id_num](i);
                    p += sizeof(double);
                  }
                break;
            }

          return p;
        };
    };

    // Integrator which chooses Euler, RK2 or RK4 depending on characteristics of the cell a particle is in
    // Currently used for research only
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

        unsigned int                    _step;
        std::map<double, Point<dim> >   _loc0, _k1, _k2, _k3;
        std::map<double, IntegrationScheme>        _scheme;
        const parallel::distributed::Triangulation<dim>   *_tria;
        const DoFHandler<dim>           *_dh;
        const Mapping<dim>              *_mapping;
        const TrilinosWrappers::MPI::BlockVector *_solution;

        virtual IntegrationScheme select_scheme(const std::vector<Point<dim> > &cell_vertices, const std::vector<Point<dim> > &cell_velocities, const double timestep)
        {
          return cell_vertices[0][0] > 0.5 ? SCHEME_RK4 : SCHEME_EULER;
        };

      public:
        HybridIntegrator(const parallel::distributed::Triangulation<dim> *new_tria, const DoFHandler<dim> *new_dh, const Mapping<dim> *new_mapping, const TrilinosWrappers::MPI::BlockVector *new_solution)
        {
          _step = 0;
          _loc0.clear();
          _k1.clear();
          _k2.clear();
          _k3.clear();
          _scheme.clear();
          _tria = new_tria;
          _dh = new_dh;
          _mapping = new_mapping;
          _solution = new_solution;
        };

        virtual bool integrate_step(std::multimap<LevelInd, T> &particles, double dt)
        {
          typename std::multimap<LevelInd, T>::iterator       it;
          Point<dim>                          loc, vel, k4;
          double                              id_num;
          LevelInd                            cur_level_ind;
          IntegrationScheme                   cur_scheme;
          //typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell;
          typename DoFHandler<dim>::active_cell_iterator found_cell;
          Functions::FEFieldFunction<dim, DoFHandler<dim>, TrilinosWrappers::MPI::BlockVector> fe_value(*_dh, *_solution, *_mapping);

          // If this is the first step, go through all the cells and determine
          // which integration scheme the particles in each cell should use
          if (_step == 0)
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
                      found_cell = typename DoFHandler<dim>::active_cell_iterator(_tria, cur_level_ind.first, cur_level_ind.second, _dh);

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
                  _scheme[it->second.get_id()] = cur_scheme;
                }
            }

          for (it=particles.begin(); it!=particles.end(); ++it)
            {
              id_num = it->second.get_id();
              loc = it->second.get_location();
              vel = it->second.get_velocity();
              switch (_scheme[id_num])
                {
                  case SCHEME_EULER:
                    if (_step == 0)
                      {
                        it->second.set_location(loc + dt*vel);
                      }
                    it->second.set_vel_check(false);
                    break;
                  case SCHEME_RK2:
                    if (_step == 0)
                      {
                        _loc0[id_num] = loc;
                        it->second.set_location(loc + 0.5*dt*vel);
                      }
                    else if (_step == 1)
                      {
                        it->second.set_location(_loc0[id_num] + dt*vel);
                      }
                    if (_step != 0) it->second.set_vel_check(false);
                    break;
                  case SCHEME_RK4:
                    if (_step == 0)
                      {
                        _loc0[id_num] = loc;
                        _k1[id_num] = dt*vel;
                        it->second.set_location(loc + 0.5*_k1[id_num]);
                      }
                    else if (_step == 1)
                      {
                        _k2[id_num] = dt*vel;
                        it->second.set_location(_loc0[id_num] + 0.5*_k2[id_num]);
                      }
                    else if (_step == 2)
                      {
                        _k3[id_num] = dt*vel;
                        it->second.set_location(_loc0[id_num] + _k3[id_num]);
                      }
                    else if (_step == 3)
                      {
                        k4 = dt*vel;
                        it->second.set_location(_loc0[id_num] + (_k1[id_num]+2*_k2[id_num]+2*_k3[id_num]+k4)/6.0);
                      }
                    break;
                  default:
                    AssertThrow(false, ExcMessage("Unknown integration scheme for hybrid integrator."));
                    break;
                }
            }

          _step = (_step+1)%4;
          if (_step == 0)
            {
              _loc0.clear();
              _k1.clear();
              _k2.clear();
              _k3.clear();
              _scheme.clear();
            }

          // Continue until we're at the last step
          return (_step != 0);
        };

        virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info)
        {
          // Add the _loc0, _k1, _k2, and _k3 data
          data_info.push_back(MPIDataInfo("loc0", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("k1", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("k2", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("k3", dim, MPI_DOUBLE, sizeof(double)));
          data_info.push_back(MPIDataInfo("scheme", 1, MPI_DOUBLE, sizeof(double)));
        };

        virtual unsigned int data_len(ParticleDataFormat format) const
        {
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                return (4*dim+1)*sizeof(double);
            }
          return 0;
        };

        virtual const char *read_data(ParticleDataFormat format, double id_num, const char *data)
        {
          const char     *p = data;
          unsigned int    i;

          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Read location data
                for (i=0; i<dim; ++i)
                  {
                    _loc0[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                // Read k1, k2 and k3
                for (i=0; i<dim; ++i)
                  {
                    _k1[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    _k2[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    _k3[id_num](i) = ((double *)p)[0];
                    p += sizeof(double);
                  }
                _scheme[id_num] = (IntegrationScheme)((double *)p)[0];
                p += sizeof(double);
                break;
            }

          return p;
        };

        virtual char *write_data(ParticleDataFormat format, double id_num, char *data)
        {
          char          *p = data;
          unsigned int  i;

          // Then write our data in the appropriate format
          switch (format)
            {
              case MPI_DATA:
              case HDF5_DATA:
                // Write location data
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _loc0[id_num](i);
                    p += sizeof(double);
                  }
                // Write k1, k2 and k3
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _k1[id_num](i);
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _k2[id_num](i);
                    p += sizeof(double);
                  }
                for (i=0; i<dim; ++i)
                  {
                    ((double *)p)[0] = _k3[id_num](i);
                    p += sizeof(double);
                  }
                // Write the integration scheme for this particle
                ((double *)p)[0] = _scheme[id_num];
                p += sizeof(double);
                break;
            }

          return p;
        };
    };
  }
}

#endif
