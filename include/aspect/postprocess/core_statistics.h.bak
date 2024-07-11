/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_core_statistics_h
#define _aspect_postprocess_core_statistics_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/boundary_temperature/dynamic_core.h>


namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that computes some statistics about the dynamic core.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class CoreStatistics : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        CoreStatistics ();
        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

        /**
         * Evaluate statistics.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

        /**
         * Export core data stored in this object. Doing this only because the boundary
         * temperature doesn't allowed to store restart data there. So we store the data needed
         * for restart here and exported to boundary temperature object if required.
         */
        const BoundaryTemperature::internal::CoreData &
        get_core_data() const;

        /**
         * Serialize the contents of this class as far as they are not read
         * from input parameter files.
         */
        template <class Archive>
        void serialize (Archive &ar, const unsigned int version);

        /**
         * Save the state of this object.
         */
        void save (std::map<std::string, std::string> &status_strings) const override;

        /**
         * Restore the state of the object.
         */
        void load (const std::map<std::string, std::string> &status_strings) override;

      private:
        /**
         * Controls whether output the total excess entropy or the individual entropy terms
         * (i.e. entropy for specific heat, radioactive heating, gravitational contribution,
         * and adiabatic contribution).
         */
        bool   excess_entropy_only;

        /**
         * Stores the core data from boundary temperature.
         */
        struct  BoundaryTemperature::internal::CoreData core_data;
    };
  }
}


#endif
