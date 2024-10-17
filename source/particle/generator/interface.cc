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

#include <aspect/particle/generator/interface.h>

#include <tuple>
#include <deal.II/grid/grid_tools.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      Interface<dim>::initialize ()
      {
        const unsigned int my_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
        random_number_generator.seed(5432+my_rank);
      }



      template <int dim>
      std::pair<Particles::internal::LevelInd,Particle<dim>>
      Interface<dim>::generate_particle(const Point<dim> &position,
                                        const types::particle_index id) const
      {
        // Try to find the cell of the given position. If the position is not
        // in the domain on the local process, throw a ExcParticlePointNotInDomain
        // exception.
        std::pair<const typename parallel::distributed::Triangulation<dim>::active_cell_iterator,
            Point<dim>> it =
              GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), position);

        // Only try to add the point if the cell it is in, is on this processor
        AssertThrow(it.first.state() == IteratorState::valid && it.first->is_locally_owned(),
                    ExcParticlePointNotInDomain());

        const Particle<dim> particle(position, it.second, id);
        const Particles::internal::LevelInd cell(it.first->level(), it.first->index());
        return std::make_pair(cell,particle);

        // Avoid warnings about missing return
        return {};
      }



      template <int dim>
      Particles::ParticleIterator<dim>
      Interface<dim>::insert_particle_at_position(const Point<dim> &position,
                                                  const types::particle_index id,
                                                  Particles::ParticleHandler<dim> &particle_handler) const
      {
        // Try to find the cell of the given position.
        const std::pair<const typename parallel::distributed::Triangulation<dim>::active_cell_iterator,
              Point<dim>> it =
                GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), position);

        if (it.first.state() != IteratorState::valid || it.first->is_locally_owned() == false)
          return particle_handler.end();

        return particle_handler.insert_particle(Particle<dim>(position, it.second, id), it.first);
      }



      template <int dim>
      std::pair<Particles::internal::LevelInd,Particle<dim>>
      Interface<dim>::generate_particle (const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell,
                                         const types::particle_index id)
      {
        // Uniform distribution on the interval [0,1]. This
        // will be used to generate random particle locations.
        std::uniform_real_distribution<double> uniform_distribution_01(0.0, 1.0);

        const BoundingBox<dim> cell_bounding_box = cell->bounding_box();

        // Generate random points in these bounds until one is within the cell
        unsigned int iteration = 0;
        const unsigned int maximum_iterations = 100;
        Point<dim> particle_position;
        while (iteration < maximum_iterations)
          {
            // First generate a random point in the bounding box...
            for (unsigned int d=0; d<dim; ++d)
              particle_position[d] = cell_bounding_box.lower_bound(d)
                                     + (uniform_distribution_01(random_number_generator) *
                                        cell_bounding_box.side_length(d));

            // ...then check whether it is actually in the cell:
            try
              {
                const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(cell, particle_position);
                if (
                  cell->reference_cell().contains_point(p_unit)
                )
                  {
                    // Add the generated particle to the set
                    const Particle<dim> new_particle(particle_position, p_unit, id);
                    const Particles::internal::LevelInd cellid(cell->level(), cell->index());
                    return std::make_pair(cellid,new_particle);
                  }
              }
            catch (typename Mapping<dim>::ExcTransformationFailed &)
              {
                // The point is not in this cell. Do nothing, just try again.
              }
            ++iteration;
          }

        // If the above algorithm has not worked (e.g. because of badly
        // deformed cells), retry generating particles
        // randomly within the reference cell. This is not generating a
        // uniform distribution in real space, but will always succeed.
        for (unsigned int d=0; d<dim; ++d)
          particle_position[d] = uniform_distribution_01(random_number_generator);

        const Point<dim> p_real = this->get_mapping().transform_unit_to_real_cell(cell,particle_position);

        // Add the generated particle to the set
        const Particle<dim> new_particle(p_real, particle_position, id);
        const Particles::internal::LevelInd cellid(cell->level(), cell->index());

        return std::make_pair(cellid, new_particle);
      }


// -------------------------------- Deal with registering models and automating
// -------------------------------- their setup and selection at run time

      namespace
      {
        std::tuple
        <aspect::internal::Plugins::UnusablePluginList,
        aspect::internal::Plugins::UnusablePluginList,
        aspect::internal::Plugins::PluginList<Interface<2>>,
        aspect::internal::Plugins::PluginList<Interface<3>>> registered_plugins;
      }



      template <int dim>
      void
      register_particle_generator (const std::string &name,
                                   const std::string &description,
                                   void (*declare_parameters_function) (ParameterHandler &),
                                   std::unique_ptr<Interface<dim>> (*factory_function) ())
      {
        std::get<dim>(registered_plugins).register_plugin (name,
                                                           description,
                                                           declare_parameters_function,
                                                           factory_function);
      }



      template <int dim>
      std::unique_ptr<Interface<dim>>
      create_particle_generator (ParameterHandler &prm)
      {
        std::string name;
        name = prm.get ("Particle generator name");

        return std::get<dim>(registered_plugins).create_plugin (name,
                                                                "Particle::Generator name");
      }



      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        // declare the entry in the parameter file
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry ("Particle generator name", "random uniform",
                           Patterns::Selection (pattern_of_names),
                           "Select one of the following models:\n\n"
                           +
                           std::get<dim>(registered_plugins).get_description_string());

        std::get<dim>(registered_plugins).declare_parameters (prm);
      }



      template <int dim>
      void
      write_plugin_graph (std::ostream &out)
      {
        std::get<dim>(registered_plugins).write_plugin_graph ("Particle generator interface",
                                                              out);
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace internal
  {
    namespace Plugins
    {
      template <>
      std::list<internal::Plugins::PluginList<Particle::Generator::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<Particle::Generator::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Generator::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<Particle::Generator::Interface<3>>::plugins = nullptr;
    }
  }

  namespace Particle
  {
    namespace Generator
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_particle_generator<dim> (const std::string &, \
                                    const std::string &, \
                                    void ( *) (ParameterHandler &), \
                                    std::unique_ptr<Interface<dim>>( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  void \
  write_plugin_graph<dim> (std::ostream &); \
  \
  template \
  std::unique_ptr<Interface<dim>> \
  create_particle_generator<dim> (ParameterHandler &prm);

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
