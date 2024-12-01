/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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

#include <aspect/particle/integrator/rk_2.h>
#include <aspect/particle/property/interface.h>
#include <aspect/particle/manager.h>
#include <aspect/geometry_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Integrator
    {
      template <int dim>
      RK2<dim>::RK2()
        :
        integrator_substep(0)
      {}



      template <int dim>
      void
      RK2<dim>::initialize ()
      {
        const auto &property_information = this->get_particle_manager(this->get_particle_manager_index()).get_property_manager().get_data_info();
        property_index_old_location = property_information.get_position_by_field_name("internal: integrator properties");
      }



      template <int dim>
      void
      RK2<dim>::local_integrate_step(const typename ParticleHandler<dim>::particle_iterator &begin_particle,
                                     const typename ParticleHandler<dim>::particle_iterator &end_particle,
                                     const std::vector<Tensor<1,dim>> &old_velocities,
                                     const std::vector<Tensor<1,dim>> &velocities,
                                     const double dt)
      {
        Assert(static_cast<unsigned int> (std::distance(begin_particle, end_particle)) == old_velocities.size(),
               ExcMessage("The particle integrator expects the old velocity vector to be of equal size "
                          "to the number of particles to advect. For some unknown reason they are different, "
                          "most likely something went wrong in the calling function."));

        if (higher_order_in_time == true && integrator_substep == 1)
          Assert(old_velocities.size() == velocities.size(),
                 ExcMessage("The particle integrator expects the velocity vector to be of equal size "
                            "to the number of particles to advect. For some unknown reason they are different, "
                            "most likely something went wrong in the calling function."));

        const auto cell = begin_particle->get_surrounding_cell();
        bool at_periodic_boundary = false;
        if (this->get_triangulation().get_periodic_face_map().empty() == false)
          for (const auto &face_index: cell->face_indices())
            if (cell->at_boundary(face_index))
              if (cell->has_periodic_neighbor(face_index))
                {
                  at_periodic_boundary = true;
                  break;
                }

        typename std::vector<Tensor<1,dim>>::const_iterator old_velocity = old_velocities.begin();
        typename std::vector<Tensor<1,dim>>::const_iterator velocity = velocities.begin();

        for (typename ParticleHandler<dim>::particle_iterator it = begin_particle;
             it != end_particle; ++it, ++velocity, ++old_velocity)
          {
            ArrayView<double> properties = it->get_properties();

            if (integrator_substep == 0)
              {
                const Tensor<1,dim> k1 = dt * (*old_velocity);
#if DEAL_II_VERSION_GTE(9, 6, 0)
                // Get a reference to the particle location, so that we can update it in-place
                Point<dim> &location = it->get_location();
#else
                Point<dim> location = it->get_location();
#endif
                Point<dim> new_location = location + 0.5 * k1;

                // Check if we crossed a periodic boundary and if necessary adjust positions
                if (at_periodic_boundary)
                  this->get_geometry_model().adjust_positions_for_periodicity(new_location,
                                                                              ArrayView<Point<dim>>(location));

                for (unsigned int i=0; i<dim; ++i)
                  {
                    properties[property_index_old_location + i] = location[i];
#if DEAL_II_VERSION_GTE(9, 6, 0)
                    location[i] = new_location[i];
#endif
                  }
#if !DEAL_II_VERSION_GTE(9, 6, 0)
                it->set_location(new_location);
#endif
              }
            else if (integrator_substep == 1)
              {
                const Tensor<1,dim> k2 = (higher_order_in_time == true)
                                         ?
                                         dt * (*old_velocity + *velocity) * 0.5
                                         :
                                         dt * (*old_velocity);

#if DEAL_II_VERSION_GTE(9, 6, 0)
                Point<dim> &location = it->get_location();
#else
                Point<dim> location = it->get_location();
#endif

                for (unsigned int i=0; i<dim; ++i)
                  location[i] = properties[property_index_old_location + i] + k2[i];

                // no need to adjust old location, because this is the last integrator step
                if (at_periodic_boundary)
                  this->get_geometry_model().adjust_positions_for_periodicity(location);

#if !DEAL_II_VERSION_GTE(9, 6, 0)
                it->set_location(location);
#endif
              }
            else
              {
                Assert(false,
                       ExcMessage("The RK2 integrator should never continue after two integration steps."));
              }
          }
      }



      template <int dim>
      bool
      RK2<dim>::new_integration_step()
      {
        integrator_substep = (integrator_substep + 1) % 2;

        // Continue until we're at the last step
        return (integrator_substep != 0);
      }



      template <int dim>
      std::array<bool, 3>
      RK2<dim>::required_solution_vectors() const
      {
        switch (integrator_substep)
          {
            case 0:
              return {{false, true, false}};
            case 1:
            {
              if (higher_order_in_time)
                return {{false, true, true}};
              else
                return {{false, true, false}};
            }
            default:
              Assert(false,
                     ExcMessage("The RK4 integrator should never continue after four integration steps."));

              return {{false, false, false}};
          }
      }



      template <int dim>
      void
      RK2<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Integrator");
        {
          prm.enter_subsection("RK2");
          {
            prm.declare_entry ("Higher order accurate in time", "true",
                               Patterns::Bool(),
                               "Whether to correctly evaluate old and current velocity "
                               "solution to reach higher-order accuracy in time. If set to "
                               "'false' only the old velocity solution is evaluated to "
                               "simulate a first order method in time. This is only "
                               "recommended for benchmark purposes.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      RK2<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Integrator");
        {
          prm.enter_subsection("RK2");
          {
            higher_order_in_time = prm.get_bool("Higher order accurate in time");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
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
      ASPECT_REGISTER_PARTICLE_INTEGRATOR(RK2,
                                          "rk2",
                                          "Second Order Runge Kutta integrator "
                                          "$y_{n+1} = y_n + \\Delta t\\, v(t_{n+1/2}, y_{n} + \\frac{1}{2} k_1)$ "
                                          "where $k_1 = \\Delta t\\, v(t_{n}, y_{n})$")
    }
  }
}
