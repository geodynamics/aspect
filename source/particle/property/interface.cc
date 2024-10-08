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

#include <aspect/particle/property/interface.h>
#include <aspect/utilities.h>

#include <aspect/particle/integrator/euler.h>
#include <aspect/particle/integrator/rk_2.h>
#include <aspect/particle/integrator/rk_4.h>


#include <aspect/boundary_composition/interface.h>

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
                                                               std::pair<std::string,unsigned int>>> &properties)
      {
        unsigned int global_component_index = 0;
        for (const auto &property : properties)
          {
            unsigned int component_per_plugin = 0;
            unsigned int field_per_plugin = 0;

            position_per_plugin.push_back(global_component_index);

            for (const auto &entry : property)
              {
                const std::string  name         = entry.first;
                const unsigned int n_components = entry.second;

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
        Assert(field_index < field_names.size(),
               ExcMessage("The number of field names (" + std::to_string(field_names.size())
                          + ") is smaller than the requested field index ("
                          + std::to_string(field_index) + ")."));

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
      void
      Interface<dim>::initialize_one_particle_property (const Point<dim> &,
                                                        std::vector<double> &) const
      {}



      DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
      template <int dim>
      void
      Interface<dim>::update_particle_properties (const ParticleUpdateInputs<dim> &inputs,
                                                  typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        Vector<double> solution;
        std::vector<Tensor<1,dim>> gradient;

        bool need_solution = false;
        bool need_gradient = false;
        unsigned int n_components = numbers::invalid_unsigned_int;

        if (inputs.solution.size() > 0)
          {
            n_components = inputs.solution[0].size();

            for (unsigned int i=0; i<n_components; ++i)
              if (get_update_flags(i) & update_values)
                {
                  need_solution = true;
                  break;
                }
          }

        if (inputs.gradients.size() > 0)
          {
            n_components = inputs.gradients[0].size();

            for (unsigned int i=0; i<n_components; ++i)
              if (get_update_flags(i) & update_gradients)
                {
                  need_gradient = true;
                  break;
                }
          }

        unsigned int i = 0;
        for (auto particle = particles.begin(); particle != particles.end(); ++particle)
          {
            if (need_solution)
              {
                solution.reinit(inputs.solution[i].size());
                for (unsigned int j=0; j<solution.size(); ++j)
                  {
                    if (get_update_flags(j) & update_values)
                      solution[j] = inputs.solution[i][j];
                    else
                      solution[j] = numbers::signaling_nan<double>();
                  }
              }

            if (need_gradient)
              {
                gradient.resize(inputs.gradients[i].size());

                for (unsigned int j=0; j<gradient.size(); ++j)
                  {
                    if (get_update_flags(j) & update_gradients)
                      gradient[j] = inputs.gradients[i][j];
                    else
                      gradient[j] = numbers::signaling_nan<Tensor<1,dim>>();
                  }
              }

            // call the deprecated version of this function
            update_particle_property(this->data_position,
                                     solution,
                                     gradient,
                                     particle);

            ++i;
          }
      }
      DEAL_II_ENABLE_EXTRA_DIAGNOSTICS



      template <int dim>
      void
      Interface<dim>::update_particle_property (const unsigned int /*data_position*/,
                                                const Vector<double> &/*solution*/,
                                                const std::vector<Tensor<1,dim>> &/*gradients*/,
                                                typename ParticleHandler<dim>::particle_iterator &/*particle*/) const
      {}



      template <int dim>
      UpdateTimeFlags
      Interface<dim>::need_update () const
      {
        return update_never;
      }



      DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
      template <int dim>
      UpdateFlags
      Interface<dim>::get_update_flags (const unsigned int /*component*/) const
      {
        // If this function is not implemented by the derived class, we use
        // the default implementation of the deprecated version of this class.
        return get_needed_update_flags();
      }
      DEAL_II_ENABLE_EXTRA_DIAGNOSTICS



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
      Interface<dim>::set_data_position (const unsigned int index)
      {
        data_position = index;
      }



      template <int dim>
      unsigned int
      Interface<dim>::get_data_position () const
      {
        return data_position;
      }



      template <int dim>
      void
      IntegratorProperties<dim>::initialize_one_particle_property(const Point<dim> &/*position*/,
                                                                  std::vector<double> &data) const
      {
        data.resize(data.size() + n_integrator_properties, 0.0);
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      IntegratorProperties<dim>::get_property_information() const
      {
        return {{"internal: integrator properties", n_integrator_properties}};
      }



      template <int dim>
      void
      IntegratorProperties<dim>::parse_parameters (ParameterHandler &prm)
      {
        std::string name;
        name = prm.get ("Integration scheme");

        if (name == "rk2")
          n_integrator_properties = Particle::Integrator::RK2<dim>::n_integrator_properties;
        else if (name == "rk4")
          n_integrator_properties = Particle::Integrator::RK4<dim>::n_integrator_properties;
        else if (name == "euler")
          n_integrator_properties = Particle::Integrator::Euler<dim>::n_integrator_properties;
        else
          AssertThrow(false,
                      ExcMessage("Unknown integrator scheme. The particle property 'Integrator properties' "
                                 "does not know how many particle properties to store for this integration scheme."));
      }



      template <int dim>
      void
      Manager<dim>::initialize ()
      {
        std::vector<std::vector<std::pair<std::string, unsigned int>>> info;

        // Get the property information of the selected plugins
        for (const auto &p : this->plugin_objects)
          {
            info.push_back(p->get_property_information());
          }

        // Initialize our property information
        property_information = ParticlePropertyInformation(info);
        unsigned int plugin_index = 0;
        for (const auto &p : this->plugin_objects)
          {
            p->set_data_position(property_information.get_position_by_plugin_index(plugin_index));
            p->initialize();
            ++plugin_index;
          }
      }



      template <int dim>
      void
      Manager<dim>::initialize_one_particle (typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        if (property_information.n_components() == 0)
          return;

        std::vector<double> particle_properties;
        particle_properties.reserve(property_information.n_components());

        for (const auto &p : this->plugin_objects)
          {
            p->initialize_one_particle_property(particle->get_location(),
                                                particle_properties);
          }

        Assert(particle_properties.size() == property_information.n_components(),
               ExcMessage("The reported numbers of particle property components do not sum up "
                          "to the number of particle properties that were initialized by "
                          "the property plugins. Check the selected property plugins for "
                          "consistency between reported size and actually set properties."));

        particle->set_properties(particle_properties);
      }



      template <int dim>
      std::vector<double>
      Manager<dim>::initialize_late_particle (const Point<dim> &particle_location,
                                              const ParticleHandler<dim> &particle_handler,
                                              const Interpolator::Interface<dim> &interpolator,
                                              const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        if (property_information.n_components() == 0)
          return {};

        std::vector<double> particle_properties;
        particle_properties.reserve(property_information.n_components());

        unsigned int property_index = 0;
        for (typename std::list<std::unique_ptr<Interface<dim>>>::const_iterator
             p = this->plugin_objects.begin(); p!=this->plugin_objects.end(); ++p, ++property_index)
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
                  (*p)->initialize_one_particle_property(particle_location,
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
                                                                                particle_location)).first;
                    }
                  else
                    found_cell = cell;

                  std::vector<std::vector<double>> interpolated_properties;

                  try
                    {
                      interpolated_properties = interpolator.properties_at_points(particle_handler,
                                                                                  std::vector<Point<dim>> (1,particle_location),
                                                                                  ComponentMask(property_information.n_components(),true),
                                                                                  found_cell);
                    }
                  // interpolators that throw exceptions usually do not result in
                  // anything good, because they result in an unwinding of the stack
                  // and, if only one processor triggers an exception, the
                  // destruction of objects often causes a deadlock or completely
                  // unrelated MPI error messages. Thus, if an exception is
                  // generated, catch it, print an error message, and abort the program.
                  catch (std::exception &exc)
                    {
                      std::cerr << std::endl << std::endl
                                << "----------------------------------------------------"
                                << std::endl;
                      std::cerr << "Exception on MPI process <"
                                << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
                                << "> while generating new particle properties: "
                                << std::endl
                                << exc.what() << std::endl
                                << "Aborting!" << std::endl
                                << "----------------------------------------------------"
                                << std::endl;

                      // terminate the program!
                      MPI_Abort (MPI_COMM_WORLD, 1);
                    }

                  for (unsigned int property_component = 0; property_component < property_information.get_components_by_plugin_index(property_index); ++property_component)
                    particle_properties.push_back(interpolated_properties[0][property_information.get_position_by_plugin_index(property_index)+property_component]);
                  break;
                }

                case aspect::Particle::Property::interpolate_respect_boundary:
                {
                  typename parallel::distributed::Triangulation<dim>::cell_iterator found_cell;

                  if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
                    {
                      found_cell = (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                                                this->get_triangulation(),
                                                                                particle_location)).first;
                    }
                  else
                    found_cell = cell;

                  const auto &manager = this->get_boundary_composition_manager();
                  const auto &fixed_boundaries = manager.get_fixed_composition_boundary_indicators();

                  // Determine if the current cell is at a Dirichlet boundary
                  bool cell_at_fixed_boundary = false;
                  unsigned int boundary_face = numbers::invalid_unsigned_int;
                  double minimum_face_distance = std::numeric_limits<double>::max();
                  for (const unsigned int f : cell->face_indices())
                    if (cell->at_boundary(f) && fixed_boundaries.count(cell->face(f)->boundary_id()) == 1)
                      {
                        const double face_center_distance = particle_location.distance_square(cell->face(f)->center(true));
                        if (face_center_distance < minimum_face_distance)
                          {
                            minimum_face_distance = face_center_distance;
                            boundary_face = f;
                            cell_at_fixed_boundary = true;
                          }
                      }

                  // If no Dirichlet boundary, interpolate
                  if (cell_at_fixed_boundary == false)
                    {
                      const std::vector<std::vector<double>> interpolated_properties = interpolator.properties_at_points(particle_handler,
                                                                                        std::vector<Point<dim>> (1,particle_location),
                                                                                        ComponentMask(property_information.n_components(),true),
                                                                                        found_cell);
                      for (unsigned int property_component = 0; property_component < property_information.get_components_by_plugin_index(property_index); ++property_component)
                        particle_properties.push_back(interpolated_properties[0][property_information.get_position_by_plugin_index(property_index)+property_component]);
                    }
                  // Otherwise use the boundary condition
                  else
                    {
                      Assert(property_information.get_components_by_plugin_index(property_index) == this->n_compositional_fields(),
                             ExcInternalError());

                      const types::boundary_id boundary_id = cell->face(boundary_face)->boundary_id();

                      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                        {
                          const double composition = manager.boundary_composition(boundary_id,particle_location,c);
                          particle_properties.push_back(composition);
                        }
                    }

                  break;
                }

                default:
                  Assert (false, ExcInternalError());
              }
          }

        Assert (particle_properties.size() == property_information.n_components(), ExcInternalError());

        return particle_properties;
      }



      template <int dim>
      void
      Manager<dim>::update_particles (ParticleUpdateInputs<dim> &inputs,
                                      typename ParticleHandler<dim>::particle_iterator_range &particles) const
      {
        unsigned int plugin_index = 0;
        for (typename std::list<std::unique_ptr<Interface<dim>>>::const_iterator
             p = this->plugin_objects.begin(); p!=this->plugin_objects.end(); ++p,++plugin_index)
          {
            (*p)->update_particle_properties(inputs,particles);
          }
      }



      template <int dim>
      UpdateTimeFlags
      Manager<dim>::need_update () const
      {
        UpdateTimeFlags update = update_never;
        for (const auto &p : this->plugin_objects)
          {
            update = std::max(update, p->need_update());
          }
        return update;
      }



      template <int dim>
      std::vector<UpdateFlags>
      Manager<dim>::get_update_flags () const
      {
        const unsigned int n_components = this->introspection().n_components;
        std::vector<UpdateFlags> update (n_components , update_default);
        for (const auto &p : this->plugin_objects)
          {
            for (unsigned int i=0; i<n_components; ++i)
              {
                update[i] |= p->get_update_flags(i);
              }
          }

        // Make sure only the flags are set that we can deal with at the moment
        for (unsigned int i=0; i<n_components; ++i)
          {
            Assert ((update[i] & ~(update_gradients | update_values)) == false,
                    ExcNotImplemented());
          }

        return update;
      }



      template <int dim>
      bool
      Manager<dim>::plugin_name_exists(const std::string &name) const
      {
        return (std::find(this->plugin_names.begin(),this->plugin_names.end(),name) != this->plugin_names.end());
      }



      template <int dim>
      bool
      Manager<dim>::check_plugin_order(const std::string &first, const std::string &second) const
      {

        AssertThrow(first != second,
                    ExcMessage("The first and second string are the same, so can not check the order."));
        AssertThrow(plugin_name_exists(first),
                    ExcMessage("Could not find a plugin with the name <" + first + ">."));
        AssertThrow(plugin_name_exists(second),
                    ExcMessage("Could not find a plugin with the name <" + second + ">."));

        return (std::find(this->plugin_names.begin(),this->plugin_names.end(),first)
                < std::find(this->plugin_names.begin(),this->plugin_names.end(),second));
      }



      template <int dim>
      unsigned int
      Manager<dim>::get_plugin_index_by_name(const std::string &name) const
      {
        const std::vector<std::string>::const_iterator plugin = std::find(this->plugin_names.begin(),
                                                                          this->plugin_names.end(),
                                                                          name);

        AssertThrow(plugin != this->plugin_names.end(),
                    ExcMessage("The particle property manager was asked for a plugin "
                               "with the name <" + name + ">, but no such plugin could "
                               "be found."));
        return std::distance(this->plugin_names.begin(),plugin);
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



      namespace
      {
        std::tuple
        <aspect::internal::Plugins::UnusablePluginList,
        aspect::internal::Plugins::UnusablePluginList,
        aspect::internal::Plugins::PluginList<Property::Interface<2>>,
        aspect::internal::Plugins::PluginList<Property::Interface<3>>> registered_plugins;
      }



      template <int dim>
      void
      Manager<dim>::declare_parameters (ParameterHandler &prm)
      {
        // finally also construct a string for Patterns::MultipleSelection that
        // contains the names of all registered particle properties
        const std::string pattern_of_names
          = std::get<dim>(registered_plugins).get_pattern_of_names ();

        prm.declare_entry("List of particle properties",
                          "",
                          Patterns::MultipleSelection(pattern_of_names),
                          "A comma separated list of particle properties that should be tracked. "
                          "By default none is selected, which means only position, velocity "
                          "and id of the particles are output. \n\n"
                          "The following properties are available:\n\n"
                          +
                          std::get<dim>(registered_plugins).get_description_string());

        // now declare the parameters of each of the registered
        // particle properties in turn
        std::get<dim>(registered_plugins).declare_parameters (prm);
      }



      template <int dim>
      void
      Manager<dim>::parse_parameters (ParameterHandler &prm)
      {
        Assert (std::get<dim>(registered_plugins).plugins != nullptr,
                ExcMessage ("No postprocessors registered!?"));

        // now also see which derived quantities we are to compute
        this->plugin_names = Utilities::split_string_list(prm.get("List of particle properties"));
        AssertThrow(Utilities::has_unique_entries(this->plugin_names),
                    ExcMessage("The list of strings for the parameter "
                               "'Particles/List of particle properties' contains entries more than once. "
                               "This is not allowed. Please check your parameter file."));

        // see if 'all' was selected (or is part of the list). if so
        // simply replace the list with one that contains all names
        if (std::find (this->plugin_names.begin(),
                       this->plugin_names.end(),
                       "all") != this->plugin_names.end())
          {
            this->plugin_names.clear();
            for (typename std::list<typename aspect::internal::Plugins::PluginList<aspect::Particle::Property::Interface<dim>>::PluginInfo>::const_iterator
                 p = std::get<dim>(registered_plugins).plugins->begin();
                 p != std::get<dim>(registered_plugins).plugins->end(); ++p)
              this->plugin_names.push_back (std::get<0>(*p));
          }

        // then go through the list, create objects and let them parse
        // their own parameters
        for (auto &plugin_name : this->plugin_names)
          {
            this->plugin_objects.emplace_back (std::get<dim>(registered_plugins)
                                               .create_plugin (plugin_name,
                                                               "Particle property plugins"));

            if (SimulatorAccess<dim> *sim = dynamic_cast<SimulatorAccess<dim>*>(this->plugin_objects.back().get()))
              sim->initialize_simulator (this->get_simulator());

            this->plugin_objects.back()->set_particle_manager_index(particle_manager_index);
            this->plugin_objects.back()->parse_parameters (prm);
          }

        // lastly store internal integrator properties:
        this->plugin_objects.emplace_back (std::make_unique<IntegratorProperties<dim>>());
        this->plugin_objects.back()->set_particle_manager_index(particle_manager_index);
        this->plugin_objects.back()->parse_parameters (prm);
        this->plugin_names.emplace_back("internal: integrator properties");
      }



      template <int dim>
      void
      Manager<dim>::set_particle_manager_index(unsigned int particle_manager_index)
      {
        // Save this value. We will tell our plugins about this, once they
        // have been created in parse_parameters().
        this->particle_manager_index = particle_manager_index;
      }



      template <int dim>
      void
      Manager<dim>::
      register_particle_property (const std::string &name,
                                  const std::string &description,
                                  void (*declare_parameters_function) (ParameterHandler &),
                                  std::unique_ptr<Property::Interface<dim>> (*factory_function) ())
      {
        std::get<dim>(registered_plugins).register_plugin (name,
                                                           description,
                                                           declare_parameters_function,
                                                           factory_function);
      }



      template <int dim>
      void
      Manager<dim>::write_plugin_graph (std::ostream &out)
      {
        std::get<dim>(registered_plugins).write_plugin_graph ("Particle property interface",
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
      std::list<internal::Plugins::PluginList<Particle::Property::Interface<2>>::PluginInfo> *
      internal::Plugins::PluginList<Particle::Property::Interface<2>>::plugins = nullptr;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Property::Interface<3>>::PluginInfo> *
      internal::Plugins::PluginList<Particle::Property::Interface<3>>::plugins = nullptr;
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

#undef INSTANTIATE
    }
  }
}
