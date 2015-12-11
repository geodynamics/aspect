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

#include <aspect/particle/property/interface.h>
#include <list>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
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
                                                        const Vector<double> &,
                                                        const std::vector<Tensor<1,dim> > &,
                                                        std::vector<double> &) const
      {}

      template <int dim>
      void
      Interface<dim>::update_one_particle_property (const unsigned int,
                                                    const Point<dim> &,
                                                    const Vector<double> &,
                                                    const std::vector<Tensor<1,dim> > &,
                                                    std::vector<double> &) const
      {}

      template <int dim>
      UpdateTimeFlags
      Interface<dim>::need_update () const
      {
        return update_never;
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
        n_property_components = 0;

        // Save the names, lengths and positions of the selected properties
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            positions.push_back(n_property_components);

            (*p)->initialize();

            const std::vector<std::pair<std::string,unsigned int> > properties = (*p)->get_property_information();
            property_component_list.insert(property_component_list.end(),properties.begin(),properties.end());

            for (unsigned int i = 0; i < properties.size(); ++i)
              {
                property_position_map.insert(std::make_pair(properties[i].first,n_property_components));
                n_property_components += properties[i].second;
              }
          }
      }

      template <int dim>
      void
      Manager<dim>::initialize_one_particle (Particle<dim> &particle,
                                             const Vector<double> &solution,
                                             const std::vector<Tensor<1,dim> > &gradients) const
      {
        std::vector<double> particle_properties (0);
        particle.set_n_property_components(n_property_components);
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            (*p)->initialize_one_particle_property(particle.get_location(),
                                                   solution,
                                                   gradients,
                                                   particle_properties);
          }
        particle.set_properties(particle_properties);
      }

      template <int dim>
      void
      Manager<dim>::update_one_particle (Particle<dim> &particle,
                                         const Vector<double> &solution,
                                         const std::vector<Tensor<1,dim> > &gradients) const
      {
        unsigned int property = 0;
        for (typename std::list<std_cxx11::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p,++property)
          {
            (*p)->update_one_particle_property(positions[property],
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
      unsigned int
      Manager<dim>::get_n_property_components () const
      {
        return n_property_components;
      }

      template <int dim>
      std::size_t
      Manager<dim>::get_particle_size () const
      {
        return (n_property_components+dim) * sizeof(double) + sizeof(types::particle_index);
      }

      template <int dim>
      const std::vector<std::pair<std::string, unsigned int> > &
      Manager<dim>::get_data_info () const
      {
        return property_component_list;
      }

      template <int dim>
      unsigned int
      Manager<dim>::get_property_component_by_name(const std::string &name) const
      {
        // see if the given name is defined
        if (property_position_map.find (name) != property_position_map.end())
          return property_position_map.find(name)->second;
        else
          {
            Assert(false,ExcMessage("The particle property manager was asked for "
                                    "the property called <" + name + ">. This "
                                    "property does not exist in this model."));
            return 0;
          }
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
          prm.enter_subsection("Tracers");
          {
            // finally also construct a string for Patterns::MultipleSelection that
            // contains the names of all registered tracer properties
            const std::string pattern_of_names
              = std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();
            prm.declare_entry("List of tracer properties",
                              "",
                              Patterns::MultipleSelection(pattern_of_names),
                              "A comma separated list of tracer properties that should be tracked. "
                              "By default none is selected, which means only position, velocity "
                              "and id of the tracers are output. \n\n"
                              "The following properties are available:\n\n"
                              +
                              std_cxx11::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        // now declare the parameters of each of the registered
        // tracer properties in turn
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
          prm.enter_subsection("Tracers");
          {
            // now also see which derived quantities we are to compute
            prop_names = Utilities::split_string_list(prm.get("List of tracer properties"));

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
