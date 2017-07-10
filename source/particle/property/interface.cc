/*
  Copyright (C) 2015 - 2017 by the authors of the ASPECT code.

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

#include <aspect/particle/property/interface.h>
#include <aspect/utilities.h>

#include <deal.II/grid/grid_tools.h>

#include <list>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ParticlePropertyInformation::ParticlePropertyInformation()
        :
        number_of_components(numbers::invalid_unsigned_int),
        number_of_fields(numbers::invalid_unsigned_int),
        number_of_plugins(numbers::invalid_unsigned_int)
      {}

      ParticlePropertyInformation::ParticlePropertyInformation(const std::vector<
                                                               std::vector<
                                                               std::pair<std::string,unsigned int> > > &properties)
      {
        unsigned int global_component_index = 0;
        for (unsigned int plugin_index = 0;
             plugin_index < properties.size(); ++plugin_index)
          {
            unsigned int component_per_plugin = 0;
            unsigned int field_per_plugin = 0;

            position_per_plugin.push_back(global_component_index);

            for (unsigned int field_index = 0;
                 field_index < properties[plugin_index].size(); ++field_index)
              {
                const std::string  name         = properties[plugin_index][field_index].first;
                const unsigned int n_components = properties[plugin_index][field_index].second;

                field_names.push_back(name);
                components_per_field.push_back(n_components);
                position_per_field.push_back(global_component_index);
                component_per_plugin += n_components;
                global_component_index += n_components;
                ++field_per_plugin;
              }

            fields_per_plugin.push_back(field_per_plugin);
            components_per_plugin.push_back(component_per_plugin);
          }

        number_of_components = global_component_index;
        number_of_fields = field_names.size();
        number_of_plugins = properties.size();
      }

      bool
      ParticlePropertyInformation::fieldname_exists(const std::string &name) const
      {
        return (std::find(field_names.begin(),field_names.end(),name) != field_names.end());
      }

      unsigned int
      ParticlePropertyInformation::get_field_index_by_name(const std::string &name) const
      {
        const std::vector<std::string>::const_iterator field = std::find(field_names.begin(),
                                                                         field_names.end(),
                                                                         name);

        AssertThrow(field != field_names.end(),
                    ExcMessage("The particle property manager was asked for a property "
                               "field with the name <" + name + ">, but no such field could "
                               "be found."));
        return std::distance(field_names.begin(),field);
      }

      std::string
      ParticlePropertyInformation::get_field_name_by_index(const unsigned int field_index) const
      {
        return field_names[field_index];
      }

      unsigned int
      ParticlePropertyInformation::get_position_by_field_name(const std::string &name) const
      {
        const unsigned int field_index = get_field_index_by_name(name);
        return position_per_field[field_index];
      }

      unsigned int
      ParticlePropertyInformation::get_components_by_field_name(const std::string &name) const
      {
        const unsigned int field_index = get_field_index_by_name(name);
        return components_per_field[field_index];
      }

      unsigned int
      ParticlePropertyInformation::get_position_by_field_index(const unsigned int field_index) const
      {
        return position_per_field[field_index];
      }

      unsigned int
      ParticlePropertyInformation::get_components_by_field_index(const unsigned int field_index) const
      {
        return components_per_field[field_index];
      }

      unsigned int
      ParticlePropertyInformation::get_position_by_plugin_index(const unsigned int plugin_index) const
      {
        return position_per_plugin[plugin_index];
      }

      unsigned int
      ParticlePropertyInformation::get_components_by_plugin_index(const unsigned int plugin_index) const
      {
        return components_per_plugin[plugin_index];
      }

      unsigned int
      ParticlePropertyInformation::get_fields_by_plugin_index(const unsigned int plugin_index) const
      {
        return fields_per_plugin[plugin_index];
      }

      unsigned int
      ParticlePropertyInformation::n_plugins() const
      {
        return number_of_plugins;
      }

      unsigned int
      ParticlePropertyInformation::n_fields() const
      {
        return number_of_fields;
      }

      unsigned int
      ParticlePropertyInformation::n_components() const
      {
        return number_of_components;
      }



      template <int dim>
      Interface<dim>::~Interface ()
      {}

      template <int dim>
      void
      Interface<dim>::initialize ()
      {}

      template <int dim>
      void
      Interface<dim>::initialize_one_particle_property (const Point<dim> &,
                                                        std::vector<double> &) const
      {}

      template <int dim>
      void
      Interface<dim>::update_one_particle_property (const unsigned int,
                                                    const Point<dim> &,
                                                    const Vector<double> &,
                                                    const std::vector<Tensor<1,dim> > &,
                                                    const ArrayView<double> &) const
      {}

      template <int dim>
      UpdateTimeFlags
      Interface<dim>::need_update () const
      {
        return update_never;
      }

      template <int dim>
      UpdateFlags
      Interface<dim>::get_needed_update_flags () const
      {
        return update_default;
      }

      template <int dim>
      InitializationModeForLateParticles
      Interface<dim>::late_initialization_mode () const
      {
        return interpolate;
      }

      template <int dim>
      void
      Interface<dim>::declare_parameters (ParameterHandler &)
      {}



      template <int dim>
      void
      Interface<dim>::parse_parameters (ParameterHandler &)
      {}

      template <int dim>
      inline
      Manager<dim>::Manager ()
      {
      }

      template <int dim>
      inline
      Manager<dim>::~Manager ()
      {
      }

      template <int dim>
      void
      Manager<dim>::initialize ()
      {
        std::vector<std::vector<std::pair<std::string, unsigned int> > > info;

        // Get the property information of the selected plugins
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            (*p)->initialize();

            info.push_back((*p)->get_property_information());
          }

        // Initialize our property information
        property_information = ParticlePropertyInformation(info);

        // Create the memory pool that will store all particle properties
        property_pool.reset(new PropertyPool(property_information.n_components()));
      }

      template <int dim>
      void
      Manager<dim>::initialize_one_particle (Particle<dim> &particle) const
      {
        if (property_information.n_components() == 0)
          return;

        std::vector<double> particle_properties;
        particle_properties.reserve(property_information.n_components());

        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            (*p)->initialize_one_particle_property(particle.get_location(),
                                                   particle_properties);
          }

        Assert(particle_properties.size() == property_information.n_components(),
               ExcMessage("The reported numbers of particle property components do not sum up "
                          "to the number of particle properties that were initialized by "
                          "the property plugins. Check the selected property plugins for "
                          "consistency between reported size and actually set properties."));

        particle.set_properties(particle_properties);
      }

      template <int dim>
      void
      Manager<dim>::initialize_late_particle (Particle<dim> &particle,
                                              const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                              const Interpolator::Interface<dim> &interpolator,
                                              const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        particle.set_property_pool(*property_pool);

        if (property_information.n_components() == 0)
          return;

        std::vector<double> particle_properties;
        particle_properties.reserve(property_information.n_components());

        unsigned int property_index = 0;
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p, ++property_index)
          {
            switch ((*p)->late_initialization_mode())
              {
                case aspect::Particle::Property::initialize_to_zero:
                {
                  for (unsigned int property_component = 0; property_component < property_information.get_components_by_plugin_index(property_index); ++property_component)
                    particle_properties.push_back(0.0);
                  break;
                }

                case aspect::Particle::Property::initialize:
                {
                  (*p)->initialize_one_particle_property(particle.get_location(),
                                                         particle_properties);
                  break;
                }

                case aspect::Particle::Property::interpolate:
                {
                  typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell;

                  if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
                    {
                      found_cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                                                this->get_triangulation(),
                                                                                particle.get_location())).first;
                    }
                  else
                    found_cell = cell;

                  const std::vector<std::vector<double> > interpolated_properties = interpolator.properties_at_points(particles,
                                                                                    std::vector<Point<dim> > (1,particle.get_location()),
                                                                                    ComponentMask(property_information.n_components(),true),
                                                                                    found_cell);
                  for (unsigned int property_component = 0; property_component < property_information.get_components_by_plugin_index(property_index); ++property_component)
                    particle_properties.push_back(interpolated_properties[0][property_information.get_position_by_plugin_index(property_index)+property_component]);
                  break;
                }

                default:
                  Assert (false, ExcInternalError());
              }
          }

        Assert (particle_properties.size() == property_information.n_components(), ExcInternalError());

        particle.set_properties(particle_properties);
      }

      template <int dim>
      void
      Manager<dim>::update_one_particle (Particle<dim> &particle,
                                         const Vector<double> &solution,
                                         const std::vector<Tensor<1,dim> > &gradients) const
      {
        unsigned int plugin_index = 0;
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p,++plugin_index)
          {
            (*p)->update_one_particle_property(property_information.get_position_by_plugin_index(plugin_index),
                                               particle.get_location(),
                                               solution,
                                               gradients,
                                               particle.get_properties());
          }
      }

      template <int dim>
      UpdateTimeFlags
      Manager<dim>::need_update () const
      {
        UpdateTimeFlags update = update_never;
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            update = std::max(update,(*p)->need_update());
          }
        return update;
      }

      template <int dim>
      UpdateFlags
      Manager<dim>::get_needed_update_flags () const
      {
        UpdateFlags update = update_default;
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            update |= (*p)->get_needed_update_flags();
          }

        return (update & (update_default | update_values | update_gradients));
      }

      template <int dim>
      unsigned int
      Manager<dim>::get_n_property_components () const
      {
        return property_information.n_components();
      }

      template <int dim>
      std::size_t
      Manager<dim>::get_particle_size () const
      {
        return (property_information.n_components()+2*dim) * sizeof(double) + sizeof(types::particle_index);
      }

      template <int dim>
      const ParticlePropertyInformation &
      Manager<dim>::get_data_info () const
      {
        return property_information;
      }

      template <int dim>
      PropertyPool &
      Manager<dim>::get_property_pool () const
      {
        return *property_pool;
      }

      template <int dim>
      unsigned int
      Manager<dim>::get_property_component_by_name(const std::string &name) const
      {
        return property_information.get_position_by_field_name(name);
      }

      namespace
      {
        std_cxx11::tuple
        <void *,
        void *,
        aspect::internal::Plugins::PluginList<Property::Interface<2> >,
        aspect::internal::Plugins::PluginList<Property::Interface<3> > > registered_plugins;
      }


      template <int dim>
      void
      Manager<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            // finally also construct a string for Patterns::MultipleSelection that
            // contains the names of all registered particle properties
            const std::string pattern_of_names
              = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();
            prm.declare_entry("List of particle properties",
                              "",
                              Patterns::MultipleSelection(pattern_of_names),
                              "A comma separated list of particle properties that should be tracked. "
                              "By default none is selected, which means only position, velocity "
                              "and id of the particles are output. \n\n"
                              "The following properties are available:\n\n"
                              +
                              std_cxx11::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        // now declare the parameters of each of the registered
        // particle properties in turn
        std_cxx11::get<dim>(registered_plugins).declare_parameters (prm);
      }


      template <int dim>
      void
      Manager<dim>::parse_parameters (ParameterHandler &prm)
      {
        Assert (std_cxx11::get<dim>(registered_plugins).plugins != 0,
                ExcMessage ("No postprocessors registered!?"));
        std::vector<std::string> prop_names;

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            // now also see which derived quantities we are to compute
            prop_names = Utilities::split_string_list(prm.get("List of particle properties"));
            AssertThrow(Utilities::has_unique_entries(prop_names),
                        ExcMessage("The list of strings for the parameter "
                                   "'Postprocess/Particles/List of particle properties' contains entries more than once. "
                                   "This is not allowed. Please check your parameter file."));

            // see if 'all' was selected (or is part of the list). if so
            // simply replace the list with one that contains all names
            if (std::find (prop_names.begin(),
                           prop_names.end(),
                           "all") != prop_names.end())
              {
                prop_names.clear();
                for (typename std::list<typename aspect::internal::Plugins::PluginList<aspect::Particle::Property::Interface<dim> >::PluginInfo>::const_iterator
                     p = std_cxx11::get<dim>(registered_plugins).plugins->begin();
                     p != std_cxx11::get<dim>(registered_plugins).plugins->end(); ++p)
                  prop_names.push_back (std_cxx11::get<0>(*p));
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        // then go through the list, create objects and let them parse
        // their own parameters
        for (unsigned int name=0; name<prop_names.size(); ++name)
          {
            aspect::Particle::Property::Interface<dim> *
            particle_property = std_cxx11::get<dim>(registered_plugins)
                                .create_plugin (prop_names[name],
                                                "Particle property plugins");

            property_list.push_back (std_cxx11::shared_ptr<Property::Interface<dim> >
                                     (particle_property));

            if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*property_list.back()))
              sim->initialize_simulator (this->get_simulator());

            property_list.back()->parse_parameters (prm);
          }
      }

      template <int dim>
      void
      Manager<dim>::
      register_particle_property (const std::string &name,
                                  const std::string &description,
                                  void (*declare_parameters_function) (ParameterHandler &),
                                  Property::Interface<dim> *(*factory_function) ())
      {
        std_cxx11::get<dim>(registered_plugins).register_plugin (name,
                                                                 description,
                                                                 declare_parameters_function,
                                                                 factory_function);
      }



      template <int dim>
      void
      Manager<dim>::write_plugin_graph (std::ostream &out)
      {
        std_cxx11::get<dim>(registered_plugins).write_plugin_graph ("Particle property interface",
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
      std::list<internal::Plugins::PluginList<Particle::Property::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Property::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Property::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Property::Interface<3> >::plugins = 0;
    }
  }

  namespace Particle
  {
    namespace Property
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  template class Manager<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}
