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



#include <aspect/particle/output/interface.h>


namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      template <int dim>
      Interface<dim>::~Interface ()
      {}

      template <int dim>
      void
      Interface<dim>::declare_parameters (ParameterHandler &)
      {}

      template <int dim>
      void
      Interface<dim>::parse_parameters (ParameterHandler &)
      {}

      template <int dim>
      void
      Interface<dim>::initialize ()
      {}

      template <int dim>
      template <class Archive>
      void Interface<dim>::serialize (Archive &ar, const unsigned int)
      {
        ar &output_file_suffix
        &file_index
        ;
      }

      template <int dim>
      void
      Interface<dim>::save (std::map<std::string, std::string> &status_strings) const
      {
        std::ostringstream os;
        aspect::oarchive oa (os);

        oa << (*this);

        status_strings["Particles output"] = os.str();
      }

      template <int dim>
      void
      Interface<dim>::load (std::map<std::string, std::string> &status_strings)
      {
        // see if something was saved
        if (status_strings.find("Particles output") != status_strings.end())
          {
            std::istringstream is (status_strings.find("Particles output")->second);
            aspect::iarchive ia (is);

            ia >> (*this);
          }
      }

      template <int dim>
      std::string
      Interface<dim>::get_file_name () const
      {
        return "particles-"
               + get_file_index()
               + "."
               + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()), 4)
               + "."
               + output_file_suffix;
      }

      template <int dim>
      std::string
      Interface<dim>::get_particle_output_location () const
      {
        return this->get_output_directory()
               + "particles/";
      }

      template <int dim>
      std::string
      Interface<dim>::get_file_index () const
      {
        return Utilities::int_to_string(file_index, 5);
      }

      template<int dim>
      void
      Interface<dim>::increment_file_index ()
      {
        file_index++;
      }

      // -------------------------------- Deal with registering models and automating
      // -------------------------------- their setup and selection at run time

      namespace
      {
        std_cxx1x::tuple
        <void *,
        void *,
        aspect::internal::Plugins::PluginList<Interface<2> >,
        aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
      }



      template <int dim>
      void
      register_particle_output (const std::string &name,
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
      create_particle_output (const std::string &name)
      {

        Interface<dim> *output
          = std_cxx1x::get<dim>(registered_plugins).create_plugin (name,
                                                                   "Particle::Output name");

        return output;
      }



      template <int dim>
      std::string
      get_names ()
      {
        return std_cxx11::get<dim>(registered_plugins).get_pattern_of_names ();
      }

      template <int dim>
      void
      declare_parameters (ParameterHandler &prm)
      {
        // declare the entry in the parameter file
        prm.enter_subsection ("Postprocess");
        {
          prm.enter_subsection ("Particles");
          {
            const std::string pattern_of_names
              = std_cxx1x::get<dim>(registered_plugins).get_pattern_of_names ();

            prm.declare_entry ("Data output format", "vtu",
                               Patterns::List (Patterns::Selection(pattern_of_names + "|none")),
                               "File format to output raw particle data in. "
                               "Multiple formats can be specified sperated by a comma (,)"
                               "If you select 'none' no output will be "
                               "written."
                               "Select one or multiple of the following models:\n\n"
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
      std::list<internal::Plugins::PluginList<Particle::Output::Interface<2> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Output::Interface<2> >::plugins = 0;
      template <>
      std::list<internal::Plugins::PluginList<Particle::Output::Interface<3> >::PluginInfo> *
      internal::Plugins::PluginList<Particle::Output::Interface<3> >::plugins = 0;
    }
  }

  namespace Particle
  {
    namespace Output
    {
#define INSTANTIATE(dim) \
  template class Interface<dim>; \
  \
  template \
  void \
  register_particle_output<dim> (const std::string &, \
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
  create_particle_output<dim> (const std::string &name);

      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}
