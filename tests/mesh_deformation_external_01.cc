/*
  Copyright (C) 2025 - 2025 by the authors of the ASPECT code.

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

#include <aspect/simulator_signals.h>
#include <aspect/simulator_access.h>

#include <iostream>

using namespace aspect;

#include <aspect/mesh_deformation/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/mesh_deformation/external_tool_interface.h>


namespace aspect
{
  namespace MeshDeformation
  {
    template <int dim>
    class TestExternalDeformation : public ExternalToolInterface<dim>
    {
      public:
        TestExternalDeformation() = default;

        void initialize() override
        {
          this->get_pcout() << "initialize()" << std::endl;
        }

        void
        update () override
        {
          this->get_pcout() << "update()" << std::endl;

          if (!this->remote_point_evaluator)
            {
              std::vector<Point<dim>> points;
              this->get_pcout() << "\tsetting points" << std::endl;

              if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator())==0)
                {
                  if constexpr (dim == 2)
                    {
                      points.emplace_back(0.3, 1.0);
                      points.emplace_back(0.5, 1.0);
                    }
                  else
                    {
                      points.emplace_back(0.3, 0.7, 1.0);
                      points.emplace_back(0.5, 0.5, 1.0);
                      points.emplace_back(0.7, 0.2, 1.0);
                      points.emplace_back(0.7, 0.91, 1.0);
                    }
                }
              else
                {
                  if constexpr (dim == 2)
                    {
                      points.emplace_back(0.72, 1.0);
                    }
                  else
                    {
                      points.emplace_back(0.9, 0.9, 1.0);
                    }
                }
              this->set_evaluation_points (points);
            }
        }

        virtual
        std::vector<Tensor<1,dim>>
        compute_updated_velocities_at_points (const std::vector<std::vector<double>> &current_solution_at_points) const override
        {
          this->get_pcout() << "compute_updated_velocities_at_points()" << std::endl;

          {
            // Copy all data to rank 0 to print to the screen.
            this->get_pcout() << "Solution at evaluation points:" << std::endl;

            std::vector<std::vector<Point<dim>>> locations_by_rank = Utilities::MPI::gather(this->get_mpi_communicator(), this->evaluation_points);
            std::vector<std::vector<std::vector<double>>> data_by_rank = Utilities::MPI::gather(this->get_mpi_communicator(), current_solution_at_points);

            const unsigned int rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
            const unsigned int size = Utilities::MPI::n_mpi_processes(this->get_mpi_communicator());

            if (rank == 0)
              {
                for (unsigned int r=0; r<size; ++r)
                  {
                    std::cout << "rank " << r << ":" << std::endl;
                    for (unsigned int index = 0; index < locations_by_rank[r].size(); ++index)
                      {
                        std::cout << locations_by_rank[r][index];
                        for (const auto &entry : data_by_rank[r][index])
                          {
                            std::cout << " " << entry;
                          }
                        std::cout << std::endl;
                      }

                  }
              }
          }

          Assert(current_solution_at_points.size() == this->evaluation_points.size(), ExcInternalError());
          std::vector<Tensor<1,dim>> velocities(current_solution_at_points.size(), Tensor<1,dim>());
          if (velocities.size()>1)
            {
              velocities[0][dim-1]=30.0;
              velocities[1][dim-1]=-3.5;
            }
          return velocities;
        }



        Tensor<1,dim>
        compute_initial_deformation_on_boundary(const types::boundary_id /*boundary_indicator*/,
                                                const Point<dim> &position) const override
        {
          const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector(position);
          Tensor<1,dim> topography_direction;
          if (gravity.norm() > 0.0)
            topography_direction = -gravity / gravity.norm();

          const double topography_amplitude = (position[0]>=0.5) ? (0.05 * (1.+std::cos(2.*numbers::PI*position[0]))) : 0.0;
          return topography_amplitude * topography_direction;
        }
    };
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(TestExternalDeformation,
                                           "external deformation",
                                           "")
  }
}
