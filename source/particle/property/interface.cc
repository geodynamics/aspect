/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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
      void
      Interface<dim>::initialize ()
      {}

      template <int dim>
      void
      Interface<dim>::initialize_particle (std::vector<double> &,
                                           const Point<dim> &,
                                           const Vector<double> &,
                                           const std::vector<Tensor<1,dim> > &)
      {}

      template <int dim>
      void
      Interface<dim>::update_particle (unsigned int &,
                                       std::vector<double> &,
                                       const Point<dim> &,
                                       const Vector<double> &,
                                       const std::vector<Tensor<1,dim> > &)
      {}

      template <int dim>
      bool
      Interface<dim>::need_update ()
      {
        return false;
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
        unsigned int old_size = 0;
        data_len = 0;

        // Save the names, lengths and positions of the selected properties
        for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            positions.push_back(data_len);

            (*p)->initialize();
            (*p)->data_names(names);
            (*p)->data_length(length);

            for (unsigned int i = old_size; i < length.size(); ++i)
              {
                property_component_map.insert(std::make_pair(names[i],data_len));
                data_len += length[i];
              }
            old_size = length.size();
          }

        // The total length is the selected properties plus the basic properties
        // position and id.
        data_len += dim + 1;
      }

      template <int dim>
      void
      Manager<dim>::initialize_particle (Particle<dim> &particle,
                                         const Vector<double> &solution,
                                         const std::vector<Tensor<1,dim> > &gradients)
      {
        std::vector<double> particle_properties (0);
        particle.set_data_len(data_len);
        for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            (*p)->initialize_particle(particle_properties,
                                      particle.get_location(),
                                      solution,
                                      gradients);
          }
        particle.set_properties(particle_properties);
      }

      template <int dim>
      void
      Manager<dim>::update_particle (Particle<dim> &particle,
                                     const Vector<double> &solution,
                                     const std::vector<Tensor<1,dim> > &gradients)
      {
        unsigned int property = 0;
        for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p,++property)
          {
            (*p)->update_particle(positions[property],
                                  particle.get_properties(),
                                  particle.get_location(),
                                  solution,
                                  gradients);
          }
      }

      template <int dim>
      bool
      Manager<dim>::need_update ()
      {
        bool update(false);
        for (typename std::list<std_cxx1x::shared_ptr<Interface<dim> > >::const_iterator
             p = property_list.begin(); p!=property_list.end(); ++p)
          {
            update = update | (*p)->need_update();
          }
        return update;
      }

      template <int dim>
      unsigned int
      Manager<dim>::get_data_len () const
      {
        return data_len;
      }

      template <int dim>
      void
      Manager<dim>::get_data_info (std::vector<std::string>  &data_names,
                                   std::vector<unsigned int> &data_length) const
      {
        data_names = names;
        data_length = length;
      }

      template <int dim>
      unsigned int
      Manager<dim>::get_property_component_by_name(const std::string &name) const
      {
        // see if the given name is defined
        if (property_component_map.find (name) != property_component_map.end())
          return property_component_map.find(name)->second;
        else
          {
            Assert(false,ExcMessage("The particle property manager was asked for "
                                    "the property called: " + name +". This property does not exist "
                                    "in this model."));
            return 0;
          }
      }

      namespace
      {
        std_cxx1x::tuple
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
              = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();
            prm.declare_entry("List of tracer properties",
                              "",
                              Patterns::MultipleSelection(pattern_of_names),
                              "A comma separated list of tracer properties that should be tracked "
                              ". By default none is selected, which means only position, velocity "
                              " and id of the tracers are outputted. \n\n"
                              "The following properties are available:\n\n"
                              +
                              std_cxx1x::get<dim>(registered_plugins).get_description_string());
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

        // now declare the parameters of each of the registered
        // tracer properties in turn
        std_cxx1x::get<dim>(registered_plugins).declare_parameters (prm);
      }


      template <int dim>
      void
      Manager<dim>::parse_parameters (ParameterHandler &prm)
      {
        Assert (std_cxx1x::get<dim>(registered_plugins).plugins != 0,
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
                     p = std_cxx1x::get<dim>(registered_plugins).plugins->begin();
                     p != std_cxx1x::get<dim>(registered_plugins).plugins->end(); ++p)
                  prop_names.push_back (std_cxx1x::get<0>(*p));
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
            particle_property = std_cxx1x::get<dim>(registered_plugins)
                                .create_plugin (prop_names[name],
                                                "Particle property plugins");

            property_list.push_back (std_cxx1x::shared_ptr<Property::Interface<dim> >
                                     (particle_property));

            if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(&*property_list.back()))
              sim->initialize (this->get_simulator());

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
        std_cxx1x::get<dim>(registered_plugins).register_plugin (name,
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
}
