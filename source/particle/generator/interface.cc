/*
  Copyright (C) 2015 - 2016 by the authors of the ASPECT code.

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

#include <deal.II/base/std_cxx1x/tuple.h>
#include <deal.II/grid/grid_tools.h>

#include <boost/lexical_cast.hpp>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      Interface<dim>::Interface()
        :
        random_number_generator(5432)
      {}

      template <int dim>
      Interface<dim>::~Interface ()
      {}

      template <int dim>
      void
      Interface<dim>::initialize ()
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
            std::pair<const typename parallel::distributed::Triangulation<dim>::active_cell_iterator,
                Point<dim> > it =
                  GridTools::find_active_cell_around_point<> (this->get_mapping(), this->get_triangulation(), position);

            // Only try to add the point if the cell it is in, is on this processor
            AssertThrow(it.first->is_locally_owned(),
                        ExcParticlePointNotInDomain());

            const Particle<dim> particle(position, it.second, id);
            const types::LevelInd cell(it.first->level(), it.first->index());
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
      std::pair<types::LevelInd,Particle<dim> >
      Interface<dim>::generate_particle (const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell,
                                         const types::particle_index id)
      {
        // Uniform distribution on the interval [0,1]. This
        // will be used to generate random particle locations.
        boost::uniform_01<double> uniform_distribution_01;

        Point<dim> max_bounds, min_bounds;
        // Get the bounds of the cell defined by the vertices
        for (unsigned int d=0; d<dim; ++d)
          {
            min_bounds[d] = std::numeric_limits<double>::max();
            max_bounds[d] = - std::numeric_limits<double>::max();
          }

        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            const Point<dim> vertex_position = cell->vertex(v);
            for (unsigned int d=0; d<dim; ++d)
              {
                min_bounds[d] = std::min(vertex_position[d], min_bounds[d]);
                max_bounds[d] = std::max(vertex_position[d], max_bounds[d]);
              }
          }

        // Generate random points in these bounds until one is within the cell
        unsigned int iteration = 0;
        const unsigned int maximum_iterations = 100;
        Point<dim> particle_position;
        while (iteration < maximum_iterations)
          {
            for (unsigned int d=0; d<dim; ++d)
              {
                particle_position[d] = uniform_distribution_01(random_number_generator) *
                                       (max_bounds[d]-min_bounds[d]) + min_bounds[d];
              }
            try
              {
                const Point<dim> p_unit = this->get_mapping().transform_real_to_unit_cell(cell, particle_position);
                if (GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                  {
                    // Add the generated particle to the set
                    const Particle<dim> new_particle(particle_position, p_unit, id);
                    const types::LevelInd cellid(cell->level(), cell->index());
                    return std::make_pair(cellid,new_particle);
                  }
              }
            catch (typename Mapping<dim>::ExcTransformationFailed &)
              {
                // The point is not in this cell. Do nothing, just try again.
              }
            iteration++;
          }
        AssertThrow (iteration < maximum_iterations,
                     ExcMessage ("Couldn't generate particle (unusual cell shape?). "
                                 "The ratio between the bounding box volume in which the particle is "
                                 "generated and the actual cell volume is approximately: " +
                                 boost::lexical_cast<std::string>(cell->measure() / (max_bounds-min_bounds).norm_square())));

        return std::make_pair(types::LevelInd(),Particle<dim>());
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
        aspect::internal::Plugins::PluginList<Interface<2> >,
        aspect::internal::Plugins::PluginList<Interface<3> > > registered_plugins;
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
          prm.enter_subsection ("Particles");
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
          prm.enter_subsection ("Particles");
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



      template <int dim>
      void
      write_plugin_graph (std::ostream &out)
      {
        std_cxx11::get<dim>(registered_plugins).write_plugin_graph ("Particle generator interface",
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
  void \
  write_plugin_graph<dim> (std::ostream &); \
  \
  template \
  Interface<dim> * \
  create_particle_generator<dim> (ParameterHandler &prm);

      ASPECT_INSTANTIATE(INSTANTIATE)
    }
  }
}

