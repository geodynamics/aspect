// ------------------------------------------------------
//
// SPDX-License-Identifier: GPL-2.0-or-later
//
// SPDX-FileCopyrightText: Copyright (C) 2015-2026 by the ASPECT authors.
//
// This file is part of ASPECT.
//
// Detailed license information governing the source code and
// contributions can be found in the folder LICENSES and in
// LICENSE and CONTRIBUTING.md at the top level directory.
//
// ------------------------------------------------------

#include <aspect/particle/population/presence_of_melt.h>
#include <aspect/particle/manager.h>

#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Particle
  {
    namespace Population
    {
      template <int dim>
      void
      PresenceOfMelt<dim>::initialize()
      {
        AssertThrow(this->introspection().compositional_name_exists("porosity"),
                    ExcMessage("The particle population model 'presence of melt' "
                               "requires a compositional field named 'porosity'."));
      }



      template <int dim>
      std::pair<unsigned int, unsigned int>
      PresenceOfMelt<dim>::particle_count_bounds(const typename DoFHandler<dim>::active_cell_iterator &cell) const
      {
        // Evaluate the porosity at quadrature points that are used for the velocity
        // when solving the Stokes system since this is where porosity is important
        // when solving the two-phase flow equations.
        const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");
        const Quadrature<dim> &quadrature_formula = this->introspection().quadratures.velocities;

        FEValues<dim> fe_values(this->get_mapping(),
                                this->get_fe(),
                                quadrature_formula,
                                update_values);
        fe_values.reinit(cell);

        std::vector<double> porosity_values(quadrature_formula.size());
        fe_values[this->introspection().extractors.compositional_fields[porosity_index]]
        .get_function_values(this->get_solution(), porosity_values);

        // If the porosity threshold is exceeded at any quadrature point, use the particle counts
        // provided by the user.
        for (const double porosity : porosity_values)
          if (porosity > porosity_threshold)
            {
              const Particle::Manager<dim> &manager = this->get_particle_manager(this->get_particle_manager_index());

              return {manager.get_min_particles_per_cell(),
                      manager.get_max_particles_per_cell()
                     };
            }
        // If the threshold is not exceeded, both bounds are zero.
        return {0, 0};
      }



      template <int dim>
      void
      PresenceOfMelt<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Population");
        {
          prm.enter_subsection("Presence of melt");
          {
            prm.declare_entry("Porosity threshold", "1e-9",
                              Patterns::Double(0),
                              "The minimum porosity required for a cell to "
                              "contain particles. If the porosity exceeds "
                              "this threshold at any quadrature point in "
                              "a cell, the minimum and maximum particle "
                              "counts specified in the Particles subsection "
                              "are used. Otherwise, both counts are zero.");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }



      template <int dim>
      void
      PresenceOfMelt<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Population");
        {
          prm.enter_subsection("Presence of melt");
          {
            porosity_threshold = prm.get_double("Porosity threshold");
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
    namespace Population
    {
      ASPECT_REGISTER_PARTICLE_POPULATION(PresenceOfMelt,
                                          "presence of melt",
                                          "A particle population model that maintains particles only "
                                          "in cells where the compositional field 'porosity' exceeds "
                                          "a prescribed threshold. Cells without melt have minimum "
                                          "and maximum particle counts of zero. These minimum and "
                                          "maximum particle counts are only enforced  if particles "
                                          "are added and removed, respectively, as part of the particle "
                                          "load balancing.")
    }
  }
}
