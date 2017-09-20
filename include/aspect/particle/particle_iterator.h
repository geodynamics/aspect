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

#ifndef _aspect_particle_particle_iterator_h
#define _aspect_particle_particle_iterator_h

#include <aspect/global.h>
#include <aspect/particle/particle_accessor.h>


namespace aspect
{
  namespace Particle
  {
    /**
     * A class that is used to iterate over particles. Together with the
     * ParticleAccessor class this is used to hide the internal implementation
     * of the particle class and the particle container.
     */
    template<int dim, int spacedim=dim>
    class ParticleIterator: public std::iterator<std::bidirectional_iterator_tag,ParticleAccessor<dim,spacedim> >
    {
      public:
        /**
         * Empty constructor. Such an object is not usable!
         */
        ParticleIterator ();

        /**
         * Constructor of the iterator. Takes a reference to the particle
         * container, and an iterator the the cell-particle pair.
         */
        ParticleIterator (const std::multimap<types::LevelInd, Particle<dim,spacedim> > &map,
                          const typename std::multimap<types::LevelInd, Particle<dim,spacedim> >::iterator &particle);

        /**
         * Dereferencing operator, returns a reference to an accessor. Usage is thus
         * like <tt>(*i).get_id ();</tt>
         */
        const ParticleAccessor<dim,spacedim> &operator * () const;

        /**
         * Dereferencing operator, non-@p const version.
         */
        ParticleAccessor<dim,spacedim> &operator * ();

        /**
         * Assignment operator.
         */
        ParticleIterator &operator = (const ParticleIterator &);

        /**
         * Dereferencing operator, returns a pointer of the particle pointed to. Usage
         * is thus like <tt>i->get_id ();</tt>
         *
         * There is a @p const and a non-@p const version.
         */
        const ParticleAccessor<dim,spacedim> *operator -> () const;

        /**
         * Dereferencing operator, non-@p const version.
         */
        ParticleAccessor<dim,spacedim> *operator -> ();

        /**
         * Compare for equality.
         */
        bool operator == (const ParticleIterator<dim,spacedim> &) const;

        /**
         * Compare for inequality.
         */
        bool operator != (const ParticleIterator<dim,spacedim> &) const;

        /**
         * Prefix <tt>++</tt> operator: <tt>++iterator</tt>. This operator advances
         * the iterator to the next element and returns a reference to
         * <tt>*this</tt>.
         */
        ParticleIterator &operator ++ ();

        /**
         * Postfix <tt>++</tt> operator: <tt>iterator++</tt>. This operator advances
         * the iterator to the next element, but returns an iterator to the element
         * previously pointed to.
         */
        ParticleIterator operator ++ (int);

        /**
         * Prefix <tt>--</tt> operator: <tt>--iterator</tt>. This operator moves
         * the iterator to the previous element and returns a reference to
         * <tt>*this</tt>.
         */
        ParticleIterator &operator -- ();

        /**
         * Postfix <tt>--</tt> operator: <tt>iterator--</tt>. This operator moves
         * the iterator to the previous element, but returns an iterator to the element
         * previously pointed to.
         */
        ParticleIterator operator -- (int);

      private:
        /**
         * The accessor to the actual particle.
         */
        ParticleAccessor<dim,spacedim> accessor;
    };
  }
}

#endif
