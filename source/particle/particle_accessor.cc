/*
  Copyright (C) 2017 by the authors of the ASPECT code.

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

#include <aspect/particle/particle_accessor.h>

namespace aspect
{
  namespace Particle
  {
    template <int dim, int spacedim>
    ParticleAccessor<dim,spacedim>::ParticleAccessor (const std::multimap<types::LevelInd, Particle<dim,spacedim> > &map,
                                                      const typename std::multimap<types::LevelInd, Particle<dim,spacedim> >::iterator &particle)
      :
      map (const_cast<std::multimap<types::LevelInd, Particle<dim,spacedim> > *> (&map)),
      particle (particle)
    {}



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::write_data (void *&data) const
    {
      Assert(particle != map->end(),
             ExcInternalError());

      particle->second.write_data(data);
    }



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::set_location (const Point<spacedim> &new_loc)
    {
      Assert(particle != map->end(),
             ExcInternalError());

      particle->second.set_location(new_loc);
    }



    template <int dim, int spacedim>
    const Point<spacedim> &
    ParticleAccessor<dim,spacedim>::get_location () const
    {
      Assert(particle != map->end(),
             ExcInternalError());

      return particle->second.get_location();
    }



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::set_reference_location (const Point<dim> &new_loc)
    {
      Assert(particle != map->end(),
             ExcInternalError());

      particle->second.set_reference_location(new_loc);
    }



    template <int dim, int spacedim>
    const Point<dim> &
    ParticleAccessor<dim,spacedim>::get_reference_location () const
    {
      Assert(particle != map->end(),
             ExcInternalError());

      return particle->second.get_reference_location();
    }



    template <int dim, int spacedim>
    types::particle_index
    ParticleAccessor<dim,spacedim>::get_id () const
    {
      Assert(particle != map->end(),
             ExcInternalError());

      return particle->second.get_id();
    }



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::set_property_pool (PropertyPool &new_property_pool)
    {
      Assert(particle != map->end(),
             ExcInternalError());

      particle->second.set_property_pool(new_property_pool);
    }



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::set_properties (const std::vector<double> &new_properties)
    {
      Assert(particle != map->end(),
             ExcInternalError());

      particle->second.set_properties(new_properties);
      return;
    }



    template <int dim, int spacedim>
    const ArrayView<const double>
    ParticleAccessor<dim,spacedim>::get_properties () const
    {
      Assert(particle != map->end(),
             ExcInternalError());

      return particle->second.get_properties();
    }



    template <int dim, int spacedim>
    typename parallel::distributed::Triangulation<dim,spacedim>::cell_iterator
    ParticleAccessor<dim,spacedim>::get_surrounding_cell (const parallel::distributed::Triangulation<dim,spacedim> &triangulation) const
    {
      Assert(particle != map->end(),
             ExcInternalError());

      const typename parallel::distributed::Triangulation<dim,spacedim>::cell_iterator cell (&triangulation,
          particle->first.first,
          particle->first.second);
      return cell;
    }



    template <int dim, int spacedim>
    const ArrayView<double>
    ParticleAccessor<dim,spacedim>::get_properties ()
    {
      Assert(particle != map->end(),
             ExcInternalError());

      return particle->second.get_properties();
    }



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::next ()
    {
      Assert (particle != map->end(),ExcInternalError());
      ++particle;
    }



    template <int dim, int spacedim>
    void
    ParticleAccessor<dim,spacedim>::prev ()
    {
      Assert (particle != map->begin(),ExcInternalError());
      --particle;
    }



    template <int dim, int spacedim>
    bool
    ParticleAccessor<dim,spacedim>::operator != (const ParticleAccessor<dim,spacedim> &other) const
    {
      return (map != other.map) || (particle != other.particle);
    }



    template <int dim, int spacedim>
    bool
    ParticleAccessor<dim,spacedim>::operator == (const ParticleAccessor<dim,spacedim> &other) const
    {
      return (map == other.map) && (particle == other.particle);
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace Particle
  {
#define INSTANTIATE(dim) \
  template class ParticleAccessor<dim>;

    ASPECT_INSTANTIATE(INSTANTIATE)
  }
}
