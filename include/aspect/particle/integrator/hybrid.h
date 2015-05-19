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

#ifndef __aspect__particle_integrator_hybrid_h
#define __aspect__particle_integrator_hybrid_h

#include <aspect/particle/integrator/interface.h>
#include <aspect/simulator_access.h>

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
      class HybridIntegrator : public Interface<dim>, SimulatorAccess<dim>
      {
        public:
        HybridIntegrator();
          virtual bool integrate_step(typename std::multimap<LevelInd, BaseParticle<dim> > &particles, const double dt);
          virtual void add_mpi_types(std::vector<MPIDataInfo> &data_info);
          virtual unsigned int data_len() const;
          virtual unsigned int read_data(const std::vector<double> &data, const unsigned int &pos, const double &id_num);
          virtual void write_data(std::vector<double> &data, const double &id_num) const;

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

          virtual IntegrationScheme select_scheme(const std::vector<Point<dim> > &cell_vertices, const std::vector<Point<dim> > &cell_velocities, const double timestep);
      };
    }
  }
}

#endif
