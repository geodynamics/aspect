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

#ifndef __aspect__particle_output_hdf5_h
#define __aspect__particle_output_hdf5_h

#include <aspect/particle/output/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/data_out_base.h>

namespace aspect
{
  namespace Particle
  {
    namespace Output
    {
      /**
       * Class that outputs particles and their properties in hdf5 format.
       *
       * @ingroup ParticleOutput
       */
      template <int dim>
      class HDF5Output : public Interface<dim>,
        public SimulatorAccess<dim>
      {
        public:
          /**
           * Constructor.
           *
           */
          HDF5Output();

          /**
           * Write data about the particles specified in the first argument
           * to a file. If possible, encode the current simulation time
           * into this file using the data provided in the last argument.
           *
           * @param[in] particles The set of particles to generate a graphical
           * representation for.
           *
           * @param [in] property_component_list A vector of the names and number
           * of components of each property. Every name entry represents the
           * name of one particle property that will be written.The number of
           * components equals one for scalar properties and dim for
           * vector properties, but any other number is valid as well
           * (e.g. number of compositional fields).
           *
           * @param[in] current_time Current time of the simulation, given as either
           * years or seconds, as selected in the input file. In other words,
           * output writers do not need to know the units in which time is
           * described.
           *
           * @return The name of the file that was written, or any other
           * information that describes what output was produced if for example
           * multiple files were created.
           */
          virtual
          std::string
          output_particle_data(const std::multimap<types::LevelInd, Particle<dim> >     &particles,
                               const std::vector<std::pair<std::string, unsigned int> > &property_component_list,
                               const double current_time);


          /**
           * Read or write the data of this object for serialization
           */
          template <class Archive>
          void serialize(Archive &ar, const unsigned int version);

          /**
           * Save the state of the object.
           */
          virtual
          void
          save (std::ostringstream &os) const;

          /**
           * Restore the state of the object.
           */
          virtual
          void
          load (std::istringstream &is);

        private:
          /**
           * Internal index of file output number.
           */
          unsigned int file_index;

          /**
           * Vector of so far created xdmf_entries (one per written output
           * file).
           */
          std::vector<XDMFEntry> xdmf_entries;
      };
    }
  }
}

#endif
