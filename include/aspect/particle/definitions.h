/*
 Copyright (C) 2011 - 2014 by the authors of the ASPECT code.

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

#ifndef __aspect__particle_definitions_h
#define __aspect__particle_definitions_h

namespace aspect
{
  namespace Particle
  {
    // TODO: in the future, upgrade multimap to ParticleMap typedef
    // with C++11 standard "using" syntax

    /// MPI tag for particle transfers
    const int           PARTICLE_XFER_TAG = 382;

    /**
     * Typedef of cell level/index pair
     */
    typedef std::pair<int, int> LevelInd;

    class MPIDataInfo
    {
      public:
        std::string     name;
        unsigned int    n_elements;

        MPIDataInfo(std::string name,
                    unsigned int num_elems)
          :
          name(name),
          n_elements(num_elems) {};
    };
  }
}

#endif
