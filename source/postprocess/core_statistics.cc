/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/postprocess/core_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/boundary_temperature/dynamic_core.h>
#include <aspect/geometry_model/spherical_shell.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    CoreStatistics<dim>::CoreStatistics()
    {
      core_data.is_initialized = false;
    }

    template <int dim>
    std::pair<std::string,std::string>
    CoreStatistics<dim>::execute (TableHandler &statistics)
    {
      // now add all of the computed heat fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;

      const BoundaryTemperature::DynamicCore<dim> &dynamic_core =
        this->get_boundary_temperature_manager().template get_matching_active_plugin<BoundaryTemperature::DynamicCore<dim>>();

      core_data = dynamic_core.get_core_data();

      // now add core mantle boundary heat flux to the statistics object
      // and create a single string that can be output to the screen
      const std::string name = "CMB heat flux out of the core (TW)";
      statistics.add_value (name, -core_data.Q/1e12);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name, 3);
      statistics.set_scientific (name, true);

      // finally have something for the screen
      screen_text.precision(3);
      screen_text << -core_data.Q/1e12 << " TW,";


      const std::string name1 = "CMB Temperature (K)";
      statistics.add_value (name1, core_data.Ti);

      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name1, 2);
      statistics.set_scientific (name1, false);

      const std::string name2 = "Inner core radius (km)";
      statistics.add_value (name2, core_data.Ri*1e-3);
      // also make sure that the other columns filled by this object
      // all show up with sufficient accuracy and in scientific notation
      statistics.set_precision (name2, 2);
      statistics.set_scientific (name2, false);

      const std::string name3 = "Light element concentration (%)";
      statistics.add_value (name3, core_data.Xi*100);
      statistics.set_precision (name3, 4);
      statistics.set_scientific (name3, false);

      if (excess_entropy_only)
        {
          const std::string name4 = "Excess entropy (W/K)";
          const double delta_E = core_data.Es*core_data.dT_dt
                                 + core_data.Er
                                 + core_data.Eh*core_data.dR_dt
                                 + core_data.El*core_data.dR_dt
                                 + core_data.Eg*core_data.dR_dt
                                 - core_data.Ek;
          statistics.add_value (name4, delta_E);
          statistics.set_precision (name4, 3);
          statistics.set_scientific (name4, true);
        }
      else
        {
          const std::string name5 = "Es (W/K)";
          statistics.add_value (name5, core_data.Es*core_data.dT_dt);
          statistics.set_precision (name5, 3);
          statistics.set_scientific (name5, true);

          const std::string name6 = "Er (W/K)";
          statistics.add_value (name6, core_data.Er);
          statistics.set_precision (name6, 3);
          statistics.set_scientific (name6, true);

          const std::string name7 = "Eh (W/K)";
          statistics.add_value (name7, core_data.Eh*core_data.dR_dt);
          statistics.set_precision (name7, 3);
          statistics.set_scientific (name7, true);

          const std::string name8 = "El (W/K)";
          statistics.add_value (name8, core_data.El*core_data.dR_dt);
          statistics.set_precision (name8, 3);
          statistics.set_scientific (name8, true);

          const std::string name9 = "Eg (W/K)";
          statistics.add_value (name9, core_data.Eg*core_data.dR_dt);
          statistics.set_precision (name9, 3);
          statistics.set_scientific (name9, true);

          const std::string name10 = "Ek (W/K)";
          statistics.add_value (name10, core_data.Ek);
          statistics.set_precision (name10, 3);
          statistics.set_scientific (name10, true);
        }

      if (dynamic_core.is_OES_used())
        {
          const std::string name11 = "Other energy source (W)";
          statistics.add_value (name11, core_data.Q_OES);
          statistics.set_precision (name11, 3);
          statistics.set_scientific (name11, true);
        }

      return std::pair<std::string, std::string> ("CMB heat flux out of the core",
                                                  screen_text.str());
    }

    template <int dim>
    void
    CoreStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Dynamic core statistics");
        {
          prm.declare_entry("Excess entropy only","false",
                            Patterns::Bool(),
                            "Output the excess entropy only instead the each entropy terms.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    CoreStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Dynamic core statistics");
        {
          excess_entropy_only = prm.get_bool("Excess entropy only");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    const BoundaryTemperature::internal::CoreData &
    CoreStatistics<dim>::get_core_data() const
    {
      return core_data;
    }

    template <int dim>
    template <class Archive>
    void CoreStatistics<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &(core_data.Ti);
      ar &(core_data.Ri);
      ar &(core_data.Xi);
      ar &(core_data.Q);
      ar &(core_data.dR_dt);
      ar &(core_data.dT_dt);
      ar &(core_data.dX_dt);
      ar &(core_data.is_initialized);
    }

    template <int dim>
    void CoreStatistics<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query the completed string below.
      std::ostringstream os;
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      status_strings["CoreStatistics"] = os.str();
    }

    template <int dim>
    void CoreStatistics<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("CoreStatistics") != status_strings.end())
        {
          std::istringstream is (status_strings.find("CoreStatistics")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(CoreStatistics,
                                  "core statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the core evolution. (Working only with dynamic core boundary temperature plugin)")
  }
}
