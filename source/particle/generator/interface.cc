/*
  Copyright (C) 2015 by the authors of the ASPECT code.

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

#include <aspect/particle/generator/interface.h>

#include <deal.II/base/std_cxx1x/tuple.h>
#include <deal.II/grid/grid_tools.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      Interface<dim>::~Interface ()
      {}

      template <int dim>
      std::pair<types::LevelInd,Particle<dim> >
      Interface<dim>::generate_particle(const Point<dim> &position,
                                        const types::particle_index id) const
      {
        // Try to find the cell of the given position. If the position is not
        // in the domain on the local process, throw a ExcParticlePointNotInDomain
        // exception.
        try
          {
            const typename parallel::distributed::Triangulation<dim>::active_cell_iterator it =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), position)).first;

            //Only try to add the point if the cell it is in, is on this processor
            AssertThrow(it->is_locally_owned(),
                        ExcParticlePointNotInDomain());

            const Particle<dim> particle(position, id);
            const types::LevelInd cell(it->level(), it->index());
            return std::make_pair(cell,particle);
          }
        catch (GridTools::ExcPointNotFound<dim> &)
          {
            AssertThrow(false,
                        ExcParticlePointNotInDomain());
          }

        // Avoid warnings about missing return
        return std::pair<types::LevelInd,Particle<dim> >();
      }

      template <int dim>
      void
      Interface<dim>::declare_parameters (ParameterHandler &)
      {}

      template <int dim>
      void
      Interface<dim>::parse_parameters (ParameterHandler &)
      {}


// -------------------------------- Deal with registering models and automating
// -------------------------------- their setup and selection at run time

      namespace
      {
        std_cxx1x::tuple
        <void *,
        void *,
        internal::Plugins::PluginList<Interface<2> >,
        internal::Plugins::PluginList<Interface<3> > > registered_plugins;
      }



      template <int dim>
      void
      register_particle_generator (const std::string &name,
                                   const std::string &description,
                                   void (*declare_parameters_function) (ParameterHandler &),
                                   Interface<dim> *(*factory_function) ())
      {
        std_cxx1x::get<dim>(registered_plugins).register_plugin (name,
                                                                 description,
                                                                 declare_parameters_function,
                                                                 factory_function);
      }


      template <int dim>
      Interface<dim> *
      create_particle_generator (ParameterHandler &prm)
      {
        std::string name;
        prm.enter_subsection ("Postprocess");
        {
          prm.enter_subsection ("Tracers");
          {
            name = prm.get ("Particle generator name");
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        return std_cxx1x::get<dim>(registered_plugins).create_plugin (name,
                                                                      "Particle::Generator name");
      }



      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        // declare the entry in the parameter file
        prm.enter_subsection ("Postprocess");
        {
          prm.enter_subsection ("Tracers");
          {
            const std::string pattern_of_names
              = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();

            prm.declare_entry ("Particle generator name", "random uniform",
                               Patterns::Selection (pattern_of_names),
                               "Select one of the following models:\n\n"
                               +
                               std_cxx1x::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();

        std_cxx1x::get<dim>(registered_plugins).declare_parameters (prm);
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
      std::list<internal::Plugins::PluginList<Particle::Generator::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Generator::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Generator::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Generator::Interface<3> >::plugins = 0;
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
                                    Interface<dim> *( *) ()); \
  \
  template  \
  void \
  declare_parameters<dim> (ParameterHandler &); \
  \
  template \
  Interface<dim> * \
  create_particle_generator<dim> (ParameterHandler &prm);

      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}

