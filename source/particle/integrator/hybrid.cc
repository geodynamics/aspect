/*
  Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/hybrid.h>
#include <deal.II/numerics/fe_field_function.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      /**
       * Integrator which chooses Euler, RK2 or RK4 depending on characteristics of the cell a particle is in.
       * Currently used for research only.
       */
        template <int dim>
          typename HybridIntegrator<dim>::IntegrationScheme
          HybridIntegrator<dim>::select_scheme(const std::vector<Point<dim> > &cell_vertices, const std::vector<Point<dim> > &cell_velocities, const double timestep)
          {
            return cell_vertices[0][0] > 0.5 ? SCHEME_RK4 : SCHEME_EULER;
          }

        template <int dim>
        HybridIntegrator<dim>::HybridIntegrator()
          {
            step = 0;
            loc0.clear();
            k1.clear();
            k2.clear();
            k3.clear();
            scheme.clear();
          }

        template <int dim>
          bool
          HybridIntegrator<dim>::integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles,
                                                const double dt)
          {
            typename std::multimap<LevelInd, BaseParticle<dim> >::iterator    it;
            const DoFHandler<dim>                            *dh = &(this->get_dof_handler());
            const Mapping<dim>                               *mapping = &(this->get_mapping());
            const parallel::distributed::Triangulation<dim>  *tria = &(this->get_triangulation());
            const LinearAlgebra::BlockVector         *solution = &(this->get_solution());
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
                          it->second.set_location(loc0[id_num] + (k1[id_num] + 2.0 * k2[id_num] + 2.0 * k3[id_num] + k4) / 6.0);
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
          }

        template <int dim>
          void
          HybridIntegrator<dim>::add_mpi_types(std::vector<MPIDataInfo> &data_info)
          {
            // Add the loc0, k1, k2, and k3 data
            data_info.push_back(MPIDataInfo("loc0", dim));
            data_info.push_back(MPIDataInfo("k1", dim));
            data_info.push_back(MPIDataInfo("k2", dim));
            data_info.push_back(MPIDataInfo("k3", dim));
            data_info.push_back(MPIDataInfo("scheme", 1));
          }

        template <int dim>
          unsigned int
          HybridIntegrator<dim>::data_len() const
          {
            return (4*dim+1);
          }

        template <int dim>
          unsigned int
          HybridIntegrator<dim>::read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num)
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
          }

        template <int dim>
          void
          HybridIntegrator<dim>::write_data(std::vector<double> &data, const double &id_num) const
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
          }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
    ASPECT_REGISTER_PARTICLE_INTEGRATOR(HybridIntegrator,
                                               "hybrid",
                                               "Integrator which chooses Euler, RK2 or RK4 depending "
                                               "on characteristics of the cell a particle is in. "
                                               "Currently used for research only.")
    }
  }
}
