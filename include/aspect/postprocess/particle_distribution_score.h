/*
  Copyright (C) 2025 - by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_particle_distribution_score_h
#define _aspect_postprocess_particle_distribution_score_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/table.h>

#include <vector>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes a score from 0 to 1 measuring
     * how clustered particles are within cells. Scores closer to 0 have more
     * even density distributions, while scores closer to 1 indicate that particles
     * are most clustered within cells.
     * @ingroup Postprocessing
     */
    template <int dim>
    class ParticleDistributionScore : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Evaluate the solution for particle clustering.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Let the postprocessor manager know about the other postprocessors
         * this one depends on. Specifically, the particles postprocessor.
         */
        std::list<std::string>
        required_other_postprocessors() const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * The `granularity` variable determines how many buckets are
         * used in the histogram which computes the density
         * distribution score.  For example, a value of 2 means
         * $2\times 2=4$ buckets in 2D.
         */
        unsigned int granularity;

        /**
         * Sorts all of the particles within the cell into a deal.II table based on their position.
         * @param cell The cell for which to compute the particle distribution.
         * @param bucket_width The size (relative to the size of the cell) of each bucket in the table.
         * @return The table with the particle information.
         */
        Table<dim,unsigned int>
        sort_particles_into_buckets(const typename Triangulation<dim>::active_cell_iterator &cell,
                                    const double bucket_width) const;
    };
  }
}


#endif
